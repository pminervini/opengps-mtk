#!/usr/bin/python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""ogps: open-source GPS daemon for the MediaTek MT6755 connsys GPS core (replaces mnld/libmnl).

usage: ogps.py [--download | --nav BRDC.rnx.gz] [--pos lat,lon,h | --cold] [--pty LINK | --nmea PATH] [--secs N]

Per chip epoch (1 Hz) it sends what stock mnld sends: FE00 next epoch, FE10 crystal model, FE01
time, one FE02 per satellite with a known orbit (predicted transmit time + Doppler), FE03 aided /
autonomous-search lists, housekeeping. Measurements (FE29) feed an EKF (pvt.py); orbits come from
IGS broadcast files or are decoded from the satellites (navmsg.py). See README.md.
"""
import argparse, calendar, collections, math, multiprocessing as mp, os, struct, sys, time
import numpy as np
import ephem, pvt, navmsg, gpscfg
from stpgps import descramble
from fwload import Chip, PLAIN, HERE, CHUNK

INIT_FE05 = ["0a0000000000000000", "3200000000", "3600000000", "2a00000000", "3400000100",
             "0f0000e803", "1500000200", "1900000300"]
OMIT = set()


def fe10(temp, c0, c1, tref=33438, step=800, valid=True):
    """FE10: crystal (DCXO) frequency model as 3 points around the chip's temperature reading
    `temp` (FE42 bytes 12-13): (temp-step, f), (temp, f), (temp+step, f), f = c0 + c1*(T - tref).
    Format: 01, then 3 x (T:u16, 0:u32, f:u16)."""
    out = b"\x01" if valid else b"\x00"           # 00 = model not trusted (mnld cold start)
    for d in (-step, 0, step):
        T = temp + d
        f = int(round(c0 + c1 * (T - tref)))
        out += struct.pack("<HIh", T & 0xFFFF, 0, f)
    return out


class NmeaPty:
    """NMEA on a pseudo-terminal for gpsd: `link` (e.g. /run/opengps/gps0) points at the pty slave.
    The slave stays open here, so the device exists whether or not gpsd is attached; when nobody
    reads, old sentences are discarded instead of blocking the daemon."""

    def __init__(self, link):
        import termios, tty
        self.termios = termios
        self.m, self.s = os.openpty()
        tty.setraw(self.s)
        os.set_blocking(self.m, False)
        os.makedirs(os.path.dirname(os.path.abspath(link)), exist_ok=True)
        if os.path.lexists(link):
            os.unlink(link)
        os.symlink(os.ttyname(self.s), link)
        os.chmod(os.ttyname(self.s), 0o666)

    def write(self, text):
        data = text.encode()
        try:
            os.write(self.m, data)
        except BlockingIOError:                    # no reader: drop what is queued, retry once
            self.termios.tcflush(self.s, self.termios.TCIFLUSH)
            try:
                os.write(self.m, data)
            except BlockingIOError:
                pass


def _ck(s):
    x = 0
    for ch in s:
        x ^= ord(ch)
    return "$%s*%02X" % (s, x)


def nmea_sentences(gpst, lat, lon, h, svs, rms):
    """GGA, RMC, GSA for one fix (geoid separation not modelled: height is ellipsoidal)."""
    u = ephem.gpst_to_unix(gpst)
    tm = time.gmtime(u)
    hms = "%02d%02d%06.3f" % (tm.tm_hour, tm.tm_min, tm.tm_sec + (u % 1))
    la, lo = abs(lat), abs(lon)
    lat_s = "%02d%07.4f,%s" % (int(la), (la - int(la)) * 60, "N" if lat >= 0 else "S")
    lon_s = "%03d%07.4f,%s" % (int(lo), (lo - int(lo)) * 60, "E" if lon >= 0 else "W")
    hdop = max(0.5, rms / 10)
    gga = _ck("GPGGA,%s,%s,%s,1,%02d,%.1f,%.1f,M,0.0,M,," % (hms, lat_s, lon_s, len(svs), hdop, h))
    rmc = _ck("GPRMC,%s,A,%s,%s,0.0,0.0,%02d%02d%02d,,,A" % (hms, lat_s, lon_s, tm.tm_mday, tm.tm_mon, tm.tm_year % 100))
    gsa = _ck("GPGSA,A,3,%s,%.1f,%.1f,%.1f" % (",".join(["%02d" % s for s in svs[:12]] + [""] * (12 - min(12, len(svs)))), hdop * 1.5, hdop, hdop))
    return [gga, rmc, gsa]


def fetch_brdc(dirn):
    """today's IGS merged broadcast ephemeris from BKG (open access); returns path or None"""
    import urllib.request
    tm = time.gmtime()
    for back in (0, 1):
        t = time.gmtime(time.time() - back * 86400)
        name = "BRDC00WRD_S_%04d%03d0000_01D_MN.rnx.gz" % (t.tm_year, t.tm_yday)
        url = "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/%04d/%03d/%s" % (t.tm_year, t.tm_yday, name)
        dst = os.path.join(dirn, name)
        try:
            if not os.path.exists(dst) or time.time() - os.path.getmtime(dst) > 3600:
                urllib.request.urlretrieve(url, dst + ".part")
                os.replace(dst + ".part", dst)
            print("broadcast ephemeris: %s" % dst, flush=True)
            return dst
        except Exception as ex:
            print("download %s failed: %s" % (url, ex), flush=True)
    return None


def _batch_job(epochs, nav, q):
    """runs in a forked process: 2+ satellite static batch solution (pvt.solve_batch2)"""
    try:
        s = pvt.solve_batch2(epochs, nav, step=6)
        q.put(s)
    except Exception as ex:        # never take the daemon down
        q.put({"error": repr(ex)})


def bringup(c, fw, cfg):
    c.send(0xFE08, cfg)
    c.wait(0xFE31)
    for i in range((len(fw) + CHUNK - 1) // CHUNK):
        idx = (i + 1).to_bytes(2, "little")
        c.send(0xFE0C, idx + fw[i * CHUNK:(i + 1) * CHUNK])
        c.wait(0xFE32, 1.0, idx)
    c.send(0xFE05, bytes.fromhex("3200000100"))
    return c.wait(0xFE42, 2.0)


SKY = {"off": None, "cand": []}      # GPS time from the satellites: confirmed offset to the system clock


def now_gpst():
    """GPS time: the system clock, corrected by the time decoded from the satellites once two
    independent subframes agree (so a wrong or unset clock does not matter)."""
    t = ephem.unix_to_gpst(time.time())
    return t + SKY["off"] if SKY["off"] is not None else t


def sky_time(ws, sv):
    """subframe words -> candidate offset (system clock vs satellites); confirmed when a second
    subframe from another satellite agrees within 1 s. Week from subframe 1 (10 bit, era 2019-2038)
    if available, else the clock's week. Returns the newly confirmed offset or None."""
    sys_t = ephem.unix_to_gpst(time.time())
    tw = navmsg.tow(ws)
    if tw >= 604800:
        return None
    week = 2048 + navmsg.bits(ws, 49, 10) if navmsg.subframe_id(ws) == 1 else int(sys_t // 604800)
    g = week * 604800 + tw
    if navmsg.subframe_id(ws) != 1 and abs(g - sys_t) > 302400:
        g += round((sys_t - g) / 604800) * 604800
    off = g - sys_t - 0.4                      # subframes reach the host ~0.35-0.5 s after their end
    SKY["cand"] = [(o, s, tt) for o, s, tt in SKY["cand"] if sys_t - tt < 60] + [(off, sv, sys_t)]
    agree = [o for o, s, tt in SKY["cand"] if s != sv and abs((o - off)) < 1.0]
    if not agree:
        return None
    off = (off + agree[-1]) / 2
    if abs(off) < 2.0:                         # the system clock is fine; pseudoranges refine it
        if SKY["off"] is not None:
            SKY["off"] = None
        return None
    if SKY["off"] is None or abs(off - SKY["off"]) > 2.0:
        SKY["off"] = off
        return off
    return None


class Daemon:
    """Chip clock model: epoch k (32-bit ms counter ms_k) happens at GPS time T_k. T is seeded
    from the system clock (NTP, ~ms), then measured every epoch from the locked channels
    (pvt.epoch_time at the known/last position). FE02 carries, per satellite, the predicted
    transmit time at T in the chip's own format (Z, hi, chips) plus the predicted Doppler."""

    def __init__(self, chip, nav, rx, log):
        self.c, self.nav, self.rx, self.log = chip, nav, np.array(rx), log
        self.started = set()
        self.epochs = 0
        self.fixes = []
        self.T = None             # GPS time of the last epoch
        self.ms = None
        self.dop_off = 0.0        # chip clock drift, Hz at L1
        self.clock_ok = False
        self.epoch_meas = []
        self.locked = 0
        self.kf = pvt.EKF(self.rx)
        self.kf_init = False
        self.cold = False             # no position: chip's autonomous search (FE03 sub 1)
        self.glo = navmsg.GloStore(nav)
        self.glo_meas = []
        self.hist = collections.deque(maxlen=240)   # cold start: (ms, locked measurements) per epoch
        self.job = None
        self.last = {}                # sv -> last FE29 (locked or pulling in)
        self.temp = 33438                 # chip temperature reading (FE42[12:14]); updated live
        self.fe10_c0, self.fe10_c1 = 416.0, 0.475
        self.fe01_drift = 0x0F8D
        self.iono_votes = collections.Counter()

    def visible(self, t):
        out = []
        for prn in self.nav["G"]:
            p = ephem.predict(self.nav, "G", prn, t, tuple(self.rx))
            if p and p["el"] > 5 and p["eph"]["health"] == 0:
                out.append((prn, p))
        return out

    def update_clock(self, ms):
        if self.T is None:
            self.T = now_gpst() - 0.08            # FE2C arrives ~80 ms after the epoch
        else:
            self.T += (ms - self.ms) / 1000.0
        self.ms = ms
        lk = [m for m in self.epoch_meas if m["state"] == 64]
        self.locked = len(lk)
        if self.cold:                  # no position yet: time from the system clock only
            return
        if getattr(self, "came_from_cold", False):   # cold fix wrong? nothing tracks -> search again
            self.no_lock = 0 if lk else getattr(self, "no_lock", 0) + 1
            if self.no_lock > 40:
                msg = "no satellite tracked for 40 s after the cold fix: back to cold search"
                print(msg, flush=True)
                self.log.write(msg + "\n")
                self.cold, self.came_from_cold, self.no_lock, self.cold_sol = True, False, 0, []
                return
            if lk:
                self.came_from_cold = len(lk) < 4      # trusted once 4+ satellites track
        if lk:
            te = pvt.epoch_time(lk, self.nav, self.rx, self.T)
            if te is not None and self.clock_ok and abs(te - self.T) > 2e-3:
                self.jumps = getattr(self, "jumps", 0) + 1
                self.log.write("clock: rejected epoch-time jump %.3f ms (%d)\n" % ((te - self.T) * 1e3, self.jumps))
                if self.jumps < 5:
                    te = None
                else:                                   # persistent: accept, restart the filter
                    self.kf = pvt.EKF(self.rx)
                    self.kf_init = False
            if te is not None and abs(te - self.T) < 0.5:
                self.jumps = 0
                dd = []
                for m in lk:
                    p = ephem.predict(self.nav, "G", m["sv"], te, tuple(self.rx))
                    if p:                                # (locked satellites without an orbit yet)
                        dd.append(m["dop"] - (-p["rate"] / ephem.C * ephem.L1))
                if not dd:
                    return
                self.log.write("clock: T correction %.3f ms, dop_off %.1f Hz, %d sv\n" % ((te - self.T) * 1e3, float(np.median(dd)), len(dd)))
                self.T, self.dop_off, self.clock_ok = te, float(np.median(dd)), True
                # closed loop on the crystal model: FE10 frequency unit ~ 1 Hz at L1, so move the
                # model offset by the measured residual; FE02 predictions keep using dop_off meanwhile
                if len(lk) >= 3:
                    step = max(-5.0, min(5.0, 0.3 * self.dop_off))   # PLLs follow only slow changes
                    self.fe10_c0 -= step                    # raising c0 raises the residual (measured: +400 -> +401 Hz)
                    self.dop_off -= step

    SWEEP = [0, 300, -300, 600, -600, 900, -900, 1200, -1200]

    def on_epoch(self, ms):
        self.update_clock(ms)
        if not self.clock_ok:          # crystal offset unknown: sweep Doppler bins until a lock
            k = (self.epochs // 3) % len(self.SWEEP)
            self.dop_off = float(self.SWEEP[k])
        self.solve()
        T, c = self.T, self.c
        nxt = ms + 1000                            # epoch counter is 32 bit: hi16, lo16
        z1, _, _ = pvt.split_tx(T + 1.0)
        c.send(0xFE00, struct.pack("<HH", nxt >> 16 & 0xFFFF, nxt & 0xFFFF))
        if "fe10" not in OMIT:
            c.send(0xFE10, fe10(self.temp, self.fe10_c0, self.fe10_c1, valid=not self.cold))
        if "fe01" not in OMIT and self.cold:      # mnld cold start: time from the clock, no epoch
            c.send(0xFE01, struct.pack("<iHHH", -3, 0x1BB4, 0, 0) + bytes([0, 1, 0]) + z1.to_bytes(2, "little") + b"\0\0")
        elif "fe01" not in OMIT:                    # [4:6] = frequency uncertainty
          c.send(0xFE01, struct.pack("<iHHH", 0, self.fe01_drift, ms >> 16 & 0xFFFF, ms & 0xFFFF)
               + bytes([(0x16 + (ms >> 16)) & 0xFF, 1, 0]) + z1.to_bytes(2, "little") + b"\0\0")
        svlist = []
        if self.cold:
            self.tail([], self.search_list())
            if self.epochs == 2:                  # mnld cold start: FE08 18 00 once, after the first lists
                c.send(0xFE08, b"\x18\x00")
            return
        for prn, p in self.visible(T):
            z, hi, chips = pvt.split_tx(T - p["range"] / ephem.C + p["clk"])
            dop = int(round(-p["rate"] / ephem.C * ephem.L1 + self.dop_off))
            first = prn not in self.started
            self.started.add(prn)
            if self.clock_ok:
                win, flags = (0x2A, 0x2E), b"\x02\x01\x01\x00"
            else:                      # (re)start acquisition each time the swept bin changes
                first = first or self.epochs % 3 == 0
                win, flags = ((0x192, 0x200), b"\x00\x00\x01\x00") if first else ((0x18C, 0x200), b"\x02\x00\x01\x00")
            c.send(0xFE02, bytes([prn, 1, 1, 0]) + struct.pack("<HBBHh", z, hi, 0, chips, dop)
                   + struct.pack("<HH", *win) + flags)
            svlist.append(bytes([prn, 1]))
        if self.glonass:
            for slot, eps in self.nav.get("R", {}).items():
                p = ephem.predict(self.nav, "R", slot, T, tuple(self.rx))
                if not p or p["el"] < 10 or p["eph"]["health"]:
                    continue
                k = p["eph"]["k"]
                z, hi, chips = pvt.split_tx_glo(T - p["range"] / ephem.C + p["clk"])
                f = 1602e6 + k * 0.5625e6
                dop = int(round(-p["rate"] / ephem.C * f + self.dop_off * f / ephem.L1))
                first = (100 + slot) not in self.started
                self.started.add(100 + slot)
                win, flags = ((0x19E, 0x100), b"\x00\x00\x01\x00") if first else ((0x69, 0x09), b"\x00\x00\x01\x00")
                c.send(0xFE02, bytes([k + 8, 0x10, 0, 0]) + struct.pack("<HBBHh", z, hi, 0, chips, dop)
                       + struct.pack("<HH", *win) + flags)
                svlist.append(bytes([k + 8, 0x10]))
        aided = {int(x[0]) for x in svlist if x[1] == 1}
        search = [x for x in self.search_list() if x not in aided][:max(0, 24 - len(svlist))] if len(svlist) < 8 else []
        self.tail(svlist, search)

    def tail(self, svlist, search=None):
        c = self.c
        c.send(0xFE03, bytes([0, len(svlist)]) + b"".join(svlist))
        if "fe03b" not in OMIT:
            # FE03 sub 1 = autonomous search list (empty when every satellite is aided)
            if search:
                c.send(0xFE03, bytes([1, 0x2E]) + b"".join(bytes([sv, 1]) for sv in search))
            else:
                c.send(0xFE03, b"\x01\x00")
        if self.epochs == 0 and "fe04" not in OMIT and not self.cold:   # mnld cold: no FE04
            c.send(0xFE04, b"\xff" * 28)

        if "fe08" not in OMIT:
            c.send(0xFE08, bytes.fromhex("060010a8e803ffffffff"))
            if self.epochs == 0:
                c.send(0xFE08, bytes.fromhex("10001000"))
        if "fe1c" not in OMIT:
            c.send(0xFE1C, b"")
        if "fe0529" not in OMIT:
            c.send(0xFE05, b"\x29\x00")
        self.epochs += 1

    # stock mnld's cold-start search order (from a capture with no aiding data)
    MNLD_ORDER = [1, 14, 29, 23, 2, 21, 11, 16, 28, 5, 20, 22, 18, 6, 19, 7, 9, 25, 8, 13, 30, 10, 3, 4]

    def search_list(self):
        """cold start as stock mnld does it: no FE02 at all; FE03 sub 1 lists up to 24 GPS PRNs
        for the chip's own blind search, locked satellites first; PRNs not in mnld's list rotate in."""
        now = time.monotonic()
        locked = [sv for sv, m in self.last.items() if m["state"] == 0x40 and now - m.get("seen", 0) < 3]
        rest = [x for x in range(1, 33) if x not in self.MNLD_ORDER]
        k = (self.epochs // 60) % len(rest)
        order = self.MNLD_ORDER[:20] + rest[k:] + rest[:k] + self.MNLD_ORDER[20:]
        return (locked + [x for x in order if x not in locked])[:24]

    def solve(self):
        ms = [m for m in self.epoch_meas if m["state"] == 64]
        self.epoch_meas = []
        glo, self.glo_meas = self.glo_meas, []
        if self.cold:
            tent = getattr(self.store, "tentative", {})
            hist_ms = [m for m in ms if m["sv"] in self.nav["G"] or m["sv"] in tent] + glo
            ms = [m for m in ms if m["sv"] in self.nav["G"]] + glo
            for m in ms:
                m["t_hint"] = now_gpst()
                m.setdefault("sys", 1)
            for m in hist_ms:
                m["t_hint"] = now_gpst()
                m.setdefault("sys", 1)
            self.hist.append((self.ms, [dict(m) for m in hist_ms]))
            if self.batch_poll():
                return
            s = None
            if len(ms) >= 3 and self.epochs % 5 == 0:
                t_rx = float(np.median([pvt.tx_from_meas(m, now_gpst()) for m in ms])) + 0.075
                s = pvt.solve_abs_alt(ms, self.nav, t_rx)
            if self.epochs % 5:
                return
            self.cold_sol = (getattr(self, "cold_sol", []) + [s["pos"]])[-3:] if s else []
            ok = len(self.cold_sol) == 3 and max(np.linalg.norm(p - self.cold_sol[-1]) for p in self.cold_sol) < 3000
            if ok:
                s["pos"] = np.mean(self.cold_sol, axis=0)
                lat, lon, h = ephem.ecef2geo(*s["pos"])
                msg = "COLD FIX %s %d sats %.6f %.6f h=%.0f rms %.1f m -> switching to aided tracking" % (
                    time.strftime("%H:%M:%S", time.gmtime()), len(s["used"]), lat, lon, h, float(np.sqrt(np.mean(s["resid"] ** 2))))
                print(msg, flush=True)
                self.log.write(msg + "\n")
                self.rx = np.array(s["pos"])
                self.cold = False
                self.clock_ok = False
                self.came_from_cold, self.no_lock = True, 0
            return
        if len(ms) < 4 or not self.clock_ok:
            return
        if not self.kf_init:                    # start the filter from a least-squares fix
            if len(ms) < 5:
                return
            s = pvt.solve(ms, self.nav, self.T, tuple(self.rx))
            if not s or float(np.sqrt(np.mean(s["resid"] ** 2))) > 50:
                return
            self.kf = pvt.EKF(s["pos"], bias_white=True)
            self.kf.ion = getattr(self, "ion", None)
            self.kf.x[6] = s["bias_s"] * ephem.C
            self.kf.P[:3, :3] = np.eye(3) * 100.0 ** 2
            self.kf.P[6, 6] = 100.0 ** 2
            self.kf_init = True
        n = self.kf.update(ms, self.nav, self.T, min_el=5.0)
        self.log.write("ekf: meas %d used %d sigma %.1f bias %.1f m\n" % (len(ms), n, float(np.sqrt(np.trace(self.kf.P[:3, :3]))), self.kf.x[6]))
        if n < 4 or np.trace(self.kf.P[:3, :3]) > 3 * 50 ** 2:
            return
        lat, lon, h = ephem.ecef2geo(*self.kf.x[:3])
        sig = float(np.sqrt(np.trace(self.kf.P[:3, :3])))
        spd = float(np.linalg.norm(self.kf.x[3:6]))
        svs = sorted(m["sv"] for m in ms)
        self.fixes.append((lat, lon, h))
        if getattr(self, "nmea", None):
            for sent in nmea_sentences(self.T, lat, lon, h, svs, sig):
                self.nmea.write(sent + "\r\n")
        msg = "FIX %s %d sats  %.6f %.6f h=%.0f  sigma %.1f m  speed %.1f m/s%s" % (
            time.strftime("%H:%M:%S", time.gmtime()), n, lat, lon, h, sig, spd, "  iono" if self.kf.ion else "")
        print(msg, flush=True)
        self.log.write(msg + "\n")

    def report_cold_fix(self, lat, lon, h):
        """a cold fix is a position fix: NMEA out, and saved as the start position for next time"""
        self.fixes.append((lat, lon, h))
        print("FIX %s cold %.6f %.6f h=%.0f" % (time.strftime("%H:%M:%S", time.gmtime()), lat, lon, h), flush=True)
        if getattr(self, "nmea", None):
            for sent in nmea_sentences(now_gpst(), lat, lon, h, [], 100.0):
                self.nmea.write(sent + "\r\n")
        try:
            import json
            json.dump(dict(lat=lat, lon=lon, h=h, t=time.time(), fe10_c0=self.fe10_c0, fe10_c1=self.fe10_c1, temp=self.temp),
                      open(self.state, "w"))
        except (OSError, AttributeError):
            pass

    def batch_poll(self):
        """start / collect the background multi-epoch solve; True if it produced the cold fix"""
        if self.job:
            proc, q = self.job
            if not q.empty():
                s = q.get()
                proc.join(1)
                self.job = None
                good = (s and "pos" in s and s["pr_rms"] < 80 and 6.3e6 < np.linalg.norm(s["pos"]) < 6.45e6
                        and len(s["svs"]) >= 3 and s["n"] >= 40 and s["span"] >= 120 and s["second"] > 2.0 * s["rms"])
                prev = getattr(self, "batch_prev", None)
                self.batch_prev = s["pos"] if good else None
                if good and (prev is None or np.linalg.norm(prev - s["pos"]) > 10e3):
                    self.log.write("batch: candidate %s (n %d span %.0f s rms %.2f second %.2f) - waiting for confirmation\n" % (
                        np.round(ephem.ecef2geo(*s["pos"]), 4), s["n"], s["span"], s["rms"], s["second"]))
                    good = False
                if good:
                    lat, lon, h = ephem.ecef2geo(*s["pos"])
                    msg = "COLD FIX (multi-epoch, %d satellites, %d measurements) %.5f %.5f h=%.0f pr rms %.1f m -> aided tracking" % (
                        len(s["svs"]), s["n"], lat, lon, h, s["pr_rms"])
                    print(msg, flush=True)
                    self.log.write(msg + "\n")
                    self.rx = np.array(s["pos"])
                    self.cold, self.clock_ok = False, False
                    self.came_from_cold, self.no_lock = False, 0      # confirmed by 2 solutions: keep it
                    self.report_cold_fix(lat, lon, h)
                    return True
                self.log.write("batch: %s\n" % (s if not s or "error" in s else "rejected pr_rms %.1f" % s["pr_rms"]))
            elif not proc.is_alive():
                self.job = None
            return False
        if self.epochs % 30 or len(self.hist) < 120:
            return False
        n = collections.Counter((m.get("sys", 1), m["sv"]) for _, ml in list(self.hist)[-180:] for m in ml)
        if sum(1 for v in n.values() if v >= 30) < 3:          # 2 satellites can be ambiguous
            return False
        q = mp.Queue()
        nav = {"G": dict(self.nav["G"]), "R": self.nav.get("R", {})}
        for sv, e in getattr(self.store, "tentative", {}).items():   # single-decoded orbits: the
            nav["G"].setdefault(sv, [e])                               # batch residuals vet them
        proc = mp.get_context("fork").Process(target=_batch_job, args=(list(self.hist)[-180:], nav, q), daemon=True)
        proc.start()
        self.job = (proc, q)
        self.log.write("batch: started with %s\n" % dict(n))
        return False

    def feed_nav(self, p):
        ws = navmsg.words(p)
        if ws[0] >> 16 == 0x8B and navmsg.subframe_id(ws) in (1, 2, 3, 4, 5) and navmsg.tow(ws) < 604800:
            off = sky_time(ws, p[2])
            if off is not None:
                msg = "TIME from the satellites: system clock %+.3f s off" % off
                print(msg, flush=True)
                self.log.write(msg + "\n")
        if ws[0] >> 16 == 0x8B:
            r = navmsg.page18(ws)
            if r:
                key = tuple(r["alpha"] + r["beta"])
                self.iono_votes[key] += 1
                k, cnt = self.iono_votes.most_common(1)[0]
                if cnt >= 2:
                    self.ion = self.kf.ion = dict(alpha=list(k[:4]), beta=list(k[4:]))
        e = self.store.feed(p, now_gpst())
        if e:
            msg = "EPH G%02d IODE %d toe %d decoded from the satellite" % (e["prn"], e["iode"], e["toe"])
            print(msg, flush=True)
            self.log.write(msg + "\n")

    def run(self, secs):
        t0 = time.monotonic()
        while not secs or time.monotonic() - t0 < secs:
            f = self.c.recv(1.0)
            if not f:
                continue
            p = f.payload if f.mtype in PLAIN else descramble(f.payload)
            self.log.write("%.3f c %04x %s\n" % (time.time(), f.mtype, p.hex()))
            if f.mtype == 0xFE42 and len(p) >= 14:
                self.temp = p[12] | p[13] << 8
            if f.mtype == 0xFE2C:
                ms = (p[0] | p[1] << 8) << 16 | (p[2] | p[3] << 8)
                self.on_epoch(ms)
                print("%s epoch ms=%6d locked=%d" % (time.strftime("%H:%M:%S"), ms, self.locked), flush=True)
            elif f.mtype == 0xFE3B and len(p) >= 14:
                e = self.glo.feed(p, now_gpst())
                if e:
                    msg = "EPH R%02d (k=%+d) tb %.0f decoded from the satellite" % (e["slot"], e["k"], e["tb"] % 86400)
                    print(msg, flush=True)
                    self.log.write(msg + "\n")
            elif f.mtype == 0xFE29 and p[2] == 0x10 and len(p) >= 33:
                m = pvt.parse_fe29(p)
                m["slot"] = self.glo.k2slot.get(p[1] - 8)
                if m["state"] == 64 and m["slot"]:
                    self.glo_meas.append(m)
            elif f.mtype == 0xFE2B and len(p) >= 46:
                self.feed_nav(p)
            elif f.mtype == 0xFE29 and p[2] == 1 and len(p) >= 33:
                m = pvt.parse_fe29(p)
                m["raw"] = p
                self.epoch_meas.append(m)
                if m["state"]:
                    m["seen"] = time.monotonic()
                    self.last[m["sv"]] = m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fw", default="/lib/firmware/mediatek/mt6755_gps.bin")
    ap.add_argument("--config", help="FE08 config as hex file (default: built by gpscfg.py)")
    ap.add_argument("--minimal", action="store_true", help="only the FE08 fields the core needs")
    ap.add_argument("--nmea", help="write NMEA (GGA/RMC/GSA) to this file or FIFO")
    ap.add_argument("--pty", help="create a pseudo-terminal for gpsd and symlink it here (e.g. /run/opengps/gps0)")
    ap.add_argument("--state", default="/var/lib/opengps/last.json", help="last position (start-up hint)")
    ap.add_argument("--download", action="store_true", help="fetch today's IGS broadcast ephemeris (BKG) if --nav is not given")
    ap.add_argument("--nav", help="RINEX broadcast nav (optional: ephemerides are also decoded from the satellites)")
    ap.add_argument("--ephcache", default="/var/lib/opengps/ephemeris.json")
    ap.add_argument("--pos", help="approximate position lat,lon,h for a warm start (default: last fix from --state; without either use --cold)")
    ap.add_argument("--secs", type=float, default=0, help="run time in seconds (0 = until stopped)")
    ap.add_argument("--fe10", help="crystal model c0,c1 for FE10 (default: learned value from --state, else 0,0.475)")
    ap.add_argument("--fe01drift", type=lambda x: int(x, 0), default=0x0F8D, help="FE01 bytes 4-5: frequency uncertainty (mnld warm ~0x0f8d, cold 0x1bb4)")
    ap.add_argument("--glonass", action="store_true", help="also track GLONASS (needs RINEX GLONASS ephemeris)")
    ap.add_argument("--debug-omit", default="", help="research: leave out host messages (comma list of fe10,fe01,init05,fe08,fe04,fe1c,fe0529,fe03b)")
    ap.add_argument("--cold", action="store_true", help="cold start: no position needed; the chip searches the sky, orbits and time come from the satellites")
    ap.add_argument("--log", help="log every chip message and solver event here (large)")
    a = ap.parse_args()
    nav = navmsg.NavStore.load(a.ephcache)
    if not a.nav and a.download:
        a.nav = fetch_brdc(os.path.dirname(a.ephcache))
    if a.nav:
        rnx = ephem.load_rinex(a.nav)
        for sv, l in rnx["G"].items():
            nav["G"].setdefault(sv, []).extend(l)
        nav["R"] = rnx["R"]
    os.makedirs(os.path.dirname(a.ephcache), exist_ok=True)
    print("ephemerides at start: %d satellites" % len(nav["G"]), flush=True)
    lat = lon = h = None
    if a.pos:
        lat, lon, h = (float(x) for x in a.pos.split(","))
    else:
        try:
            import json
            st = json.load(open(a.state))
            lat, lon, h = st["lat"], st["lon"], st["h"]
            print("start position from %s: %.4f %.4f" % (a.state, lat, lon), flush=True)
        except (OSError, ValueError, KeyError):
            pass
    if lat is None:
        if not a.cold:
            print("no position (--pos or a last fix in --state): doing a cold start", flush=True)
        a.cold = True
        lat, lon, h = 0.0, 0.0, 0.0
    fw = open(a.fw, "rb").read()
    cfg = bytes.fromhex("".join(l.split("#")[0] for l in open(a.config))) if a.config else gpscfg.build(a.minimal)
    c = Chip()
    bringup(c, fw, cfg)
    OMIT.update(x for x in a.debug_omit.split(",") if x)
    if "init05" not in OMIT:
        for hx in INIT_FE05:
            if a.cold and hx == "1900000300":      # mnld's cold init has no 0x19 word
                continue
            c.send(0xFE05, bytes.fromhex(hx))
    print("GPS core up; driving it without mnld", flush=True)
    d = Daemon(c, nav, ephem.geo2ecef(lat, lon, h), open(a.log or os.devnull, "w", buffering=1))
    d.store = navmsg.NavStore(nav, a.ephcache)
    d.nmea = NmeaPty(a.pty) if a.pty else open(a.nmea, "w", buffering=1) if a.nmea else None
    d.state = a.state
    d.glonass = a.glonass
    d.cold = a.cold
    c0c1 = a.fe10
    if not c0c1:
        try:
            import json
            st = json.load(open(a.state))
            c0c1 = "%f,%f" % (st["fe10_c0"], st["fe10_c1"])
        except (OSError, ValueError, KeyError):
            c0c1 = "0,0.475"
    d.fe10_c0, d.fe10_c1 = (float(x) for x in c0c1.split(","))
    print("crystal model: c0 %.1f Hz, c1 %.3f Hz/LSB" % (d.fe10_c0, d.fe10_c1), flush=True)
    d.fe01_drift = a.fe01drift
    try:
        d.run(a.secs)
    finally:
        os.close(c.fd)
    if d.fixes:
        arr = np.array(d.fixes)
        try:
            import json
            json.dump(dict(lat=float(arr[-1][0]), lon=float(arr[-1][1]), h=float(arr[-1][2]), t=time.time(),
                           fe10_c0=d.fe10_c0 + d.fe10_c1 * (33438 - 33438), fe10_c1=d.fe10_c1, temp=d.temp), open(a.state, "w"))
        except OSError:
            pass
        print("fixes: %d, median %.6f %.6f h=%.0f" % (len(arr), *np.median(arr, axis=0)))
    else:
        print("no fix")


if __name__ == "__main__":
    main()
