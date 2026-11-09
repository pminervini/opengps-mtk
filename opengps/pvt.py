# SPDX-License-Identifier: GPL-2.0-or-later
"""Position/velocity/time from FE29 measurements + broadcast ephemeris (GPS L1 C/A).

FE29 (descrambled, 57 bytes) fields used here:
  [0] channel  [1] sv  [2] system (01 GPS, 10 GLONASS: sv = frequency channel + 8)  [3] C/N0 dB-Hz
  [6] tracking state flags   [23:25] Z-count (1.5 s, low 16 bits of TOW/1.5)  [25] sub-tick (62.5 ms)
  [27:29] transmit time, integer chips mod 2^16   [29:31] fraction of a chip / 65536
  [31:33] Doppler, Hz (int16)
The receive time is the chip's own clock (Z*1.5 + sub*0.0625 s), so a common clock bias is solved.
"""
import math
import numpy as np
import ephem

CHIP = 1 / 1.023e6
WRAP = 65536 * CHIP                      # 64.06 ms
LD = np.longdouble                       # transmit times: absolute GPS seconds (~1.5e9) carry only
                                         # ~2.4e-7 s (72 m) resolution in float64 - keep them in
                                         # extended precision; receive-time rounding is common to all
                                         # satellites of an epoch and is absorbed by the clock bias


def parse_fe29(p):
    return dict(ch=p[0], sv=p[1], sys=p[2], cn0=p[3], state=p[6], flags=p[5], z=p[23] | p[24] << 8,
                sub=p[25], tx_chips=(p[27] | p[28] << 8) + (p[29] | p[30] << 8) / 65536.0,
                dop=int.from_bytes(p[31:33], "little", signed=True))


def full_rx_time(z, sub, t_hint):
    """chip tag -> GPS seconds (gpst), using t_hint (gpst, +-12 h) to resolve the Z-count wrap."""
    week0 = t_hint // 604800 * 604800
    zf = round((t_hint - week0) / 1.5)
    zf -= ((zf - z + 32768) % 65536) - 32768
    return week0 + zf * 1.5 + sub * 0.0625


def tropo(el_deg):
    return 2.47 / (math.sin(math.radians(max(el_deg, 3))) + 0.0121)


def z_start(t_rx, sub):
    """start of the Z-count interval the chip references its transmit-time count to"""
    return t_rx - sub * 0.0625


def tx_from_meas(m, t_hint):
    if "tx_override" in m:
        return m["tx_override"]
    return _tx_from_meas(m, t_hint)


def _tx_from_meas(m, t_hint):
    """satellite (SV-clock) transmit time from FE29 bytes 23..30:
    Z-count * 1.5 s + (hi * 65536 + chips + frac) chips. hi ([25], 0..23) is the high byte of
    the chip count since the start of the Z interval (1.5 s = 1534500 chips)."""
    return LD(full_rx_time(m["z"], 0, t_hint)) + LD(m["sub"] * 65536 + m["tx_chips"]) * LD(CHIP)


def split_tx(tx):
    """GPS transmit time -> (Z low 16 bits, hi byte, chips u16) as used by FE02/FE29"""
    tow = tx % 604800
    z = int(tow // 1.5)
    n = int(round((tow - z * 1.5) * 1.023e6))
    return z & 0xFFFF, n >> 16, n & 0xFFFF


def epoch_time(meas, nav, x, t_hint):
    """GPS time of the receive epoch from locked channels at known position x (median)."""
    v = []
    for m in meas:
        tx = tx_from_meas(m, t_hint)
        e = ephem.pick(nav["G"].get(m["sv"], []), tx, "toe_g")
        if not e:
            continue
        p = ephem.predict(nav, "G", m["sv"], tx + 0.075, tuple(x))
        if not p:
            continue
        v.append(tx + p["range"] / ephem.C - p["clk"])
    return float(np.median(v)) if v else None


def _wrap(x, p):
    return (x + p / 2) % p - p / 2


def solve(meas, nav, t_rx, x0, min_el=5.0, iters=10):
    """t_rx: one receive time for all, or None to use each measurement's own m["t_rx"]."""
    return _solve([(m, t_rx if t_rx is not None else m["t_rx"]) for m in meas], nav, x0, min_el, iters)


def _solve(mt, nav, x0, min_el=5.0, iters=10):
    """meas: parse_fe29 dicts (GPS, same epoch), t_rx: chip tag time (gpst), x0: approximate ECEF
    position (within ~100 km). Uses only the code phase (mod 1 ms) and resolves whole ms from the
    predicted ranges, so bit/frame ambiguities do not matter.
    Returns dict(pos, bias_s, resid, used) or None."""
    x0 = np.array(x0, float)
    pre = []
    for m, t_rx in mt:
        e = ephem.pick(nav["G"].get(m["sv"], []), t_rx, "toe_g")
        if not e:
            continue
        p = ephem.predict(nav, "G", m["sv"], t_rx, tuple(x0))
        if not p or p["el"] < min_el:
            continue
        tpred = t_rx - p["range"] / ephem.C + p["clk"]
        d = _wrap(tx_from_meas(m, t_rx) - tpred, 1e-3)
        pre.append((m, e, tpred, d, p["el"], t_rx))
    if len(pre) < 4:
        return None
    ang = np.array([q[3] for q in pre]) / 1e-3 * 2 * np.pi
    common = np.angle(np.mean(np.exp(1j * ang))) / (2 * np.pi) * 1e-3
    obs = []
    for m, e, tpred, d, el, t_rx in pre:
        txm = tx_from_meas(m, t_rx)                              # measured, mod 1 ms
        k = round(float((tpred + common + _wrap(d - common, 1e-3)) - txm) / 1e-3)
        tx = txm + LD(k) * LD(1e-3)                              # sv time, measured precision
        obs.append((m, e, tx, el, t_rx))
    x = x0.copy()
    b = 0.0
    for it in range(iters):
        H, r, used = [], [], []
        for m, e, tx, el, t_rx in obs:
            pos, vel, clk, _ = ephem.gps_sat(e, tx)
            tx_true = tx - clk
            pos, vel, clk, _ = ephem.gps_sat(e, tx_true)
            th = ephem.OMEGA_E * (t_rx - tx_true)
            sp = np.array([pos[0] * math.cos(th) + pos[1] * math.sin(th),
                           -pos[0] * math.sin(th) + pos[1] * math.cos(th), pos[2]])
            rho = float(LD(ephem.C) * (LD(t_rx) - tx)) + ephem.C * clk - tropo(el)
            d = sp - x
            rng = np.linalg.norm(d)
            H.append(np.r_[-d / rng, 1.0])
            r.append(rho - rng - b)
            used.append((m["sv"], el))
        H, r = np.array(H), np.array(r)
        dx = np.linalg.lstsq(H, r, rcond=None)[0]
        x += dx[:3]
        b += dx[3]
        if np.linalg.norm(dx[:3]) < 1e-3:
            break
    return dict(pos=x, bias_s=b / ephem.C, resid=r - H @ dx, used=used)


LAM = ephem.C / ephem.L1


class EKF:
    """Single-point GPS EKF (RTKLIB 'single' + filter): state x,y,z,vx,vy,vz,clock bias b (m),
    drift d (m/s). Measurements: pseudorange (iono Klobuchar, tropo, elevation/CN0 weighting) and
    Doppler-derived range rate. Static-friendly process noise; outliers rejected by innovation."""

    def __init__(self, x0, q_pos=0.05, q_clk=100.0, bias_white=False):
        """bias_white: the epoch time is re-anchored to the satellites every epoch (ogps), so the
        clock bias is re-estimated each epoch rather than propagated."""
        self.bias_white = bias_white
        self.x = np.zeros(8)
        self.x[:3] = x0
        self.P = np.diag([1e4] * 3 + [25.0] * 3 + [1e12, 1e4])
        self.q_pos, self.q_clk = q_pos, q_clk
        self.t = None
        self.ion = None

    def predict(self, t):
        if self.t is None:
            self.t = t
            return
        dt = t - self.t
        self.t = t
        F = np.eye(8)
        F[0:3, 3:6] = np.eye(3) * dt
        F[6, 7] = dt
        Q = np.zeros((8, 8))
        Q[3:6, 3:6] = np.eye(3) * self.q_pos * dt
        Q[0:3, 0:3] = np.eye(3) * self.q_pos * dt ** 3 / 3
        Q[6, 6] = 1e10 if self.bias_white else self.q_clk * dt
        Q[7, 7] = self.q_clk * dt
        self.x = F @ self.x
        self.P = F @ self.P @ F.T + Q

    def update(self, meas, nav, t_rx, min_el=10.0, gate=5.0):
        """meas: parse_fe29 dicts (locked GPS) of one epoch; t_rx: epoch GPS time (chip clock model).
        Returns number of measurements used."""
        self.predict(t_rx)
        x = self.x[:3]
        lat, lon, _ = ephem.ecef2geo(*x)
        rows, z, h, R = [], [], [], []
        for m in meas:
            e = ephem.pick(nav["G"].get(m["sv"], []), t_rx, "toe_g")
            if not e:
                continue
            p = ephem.predict(nav, "G", m["sv"], t_rx, tuple(x))
            if not p or p["el"] < min_el:
                continue
            tpred = t_rx - p["range"] / ephem.C + p["clk"]
            txm = tx_from_meas(m, t_rx)
            k = round(float(tpred + _wrap(float(txm - tpred) - self.x[6] / ephem.C, 1e-3) + self.x[6] / ephem.C - txm) / 1e-3)
            tx = txm + LD(k) * LD(1e-3)
            pos, vel, clk, _ = ephem.gps_sat(e, tx)
            pos, vel, clk, drift = ephem.gps_sat(e, tx - clk)
            th = ephem.OMEGA_E * (t_rx - tx + clk)
            sp = np.array([pos[0] * math.cos(th) + pos[1] * math.sin(th), -pos[0] * math.sin(th) + pos[1] * math.cos(th), pos[2]])
            sv = np.array(vel)
            d = sp - x
            rng = np.linalg.norm(d)
            u = d / rng
            iono = klob(self.ion, lat, lon, p["az"], p["el"], t_rx) if self.ion else 0.0
            rho = float(LD(ephem.C) * (LD(t_rx) - tx)) + ephem.C * clk - tropo(p["el"]) - iono
            s2 = (0.3 ** 2 + 0.3 ** 2 / math.sin(math.radians(p["el"])) ** 2) * 10 ** ((40 - min(m["cn0"], 40)) / 10) * 25
            rows.append(np.r_[-u, 0, 0, 0, 1, 0]); z.append(rho); h.append(rng + self.x[6]); R.append(s2)
            # Doppler -> pseudorange rate
            rr = -m["dop"] * LAM + ephem.C * drift
            rows.append(np.r_[0, 0, 0, -u, 0, 1]); z.append(rr); h.append(u @ (sv - self.x[3:6]) + self.x[7]); R.append(0.5 ** 2 * 10 ** ((40 - min(m["cn0"], 40)) / 10))
        if len(z) < 8:
            return 0
        H, z, h, R = np.array(rows), np.array(z), np.array(h), np.diag(R)
        if self.t is not None and np.trace(self.P[:3, :3]) > 1e6:        # first fix: snap clock bias
            v = z[0::2] - h[0::2]
            self.x[6] += np.median(v)
            h[0::2] += np.median(v)
        # innovation gating
        S = H @ self.P @ H.T + R
        y = z - h
        ok = np.abs(y) < gate * np.sqrt(np.abs(np.diag(S))) + 50
        H, y, R = H[ok], y[ok], R[np.ix_(ok, ok)]
        if len(y) < 8:
            return 0
        S = H @ self.P @ H.T + R
        K = self.P @ H.T @ np.linalg.inv(S)
        self.x = self.x + K @ y
        I_KH = np.eye(8) - K @ H
        self.P = I_KH @ self.P @ I_KH.T + K @ R @ K.T          # Joseph form
        self.P = (self.P + self.P.T) / 2
        return int(ok.sum() // 2)


def klob(ion, lat, lon, az, el, t):
    import navmsg
    return navmsg.klobuchar(ion, lat, lon, az, el, t % 86400)


GLO_CHIP = 1 / 0.511e6


def glo_tod(gpst):
    """GPS time -> GLONASS time of day (UTC(SU) + 3 h), seconds"""
    return (ephem.gpst_to_unix(gpst) + 3 * 3600) % 86400


def split_tx_glo(tx_gpst):
    """GLONASS transmit time (as GPS seconds) -> (TOD s low 16 bits, hi, chips u16) for FE02/FE29"""
    tod = glo_tod(tx_gpst)
    s = int(tod)
    n = int(round((tod - s) * 0.511e6))
    return s & 0xFFFF, n >> 16, n & 0xFFFF


def solve_abs(meas, nav, t_rx, iters=12):
    """Cold-start fix: no prior position. Uses full transmit times (needs locked channels with
    TOW sync), starts from the Earth's centre. Returns dict(pos, bias_s, resid, used) or None."""
    obs = []
    for m in meas:
        e = ephem.pick(nav["G"].get(m["sv"], []), t_rx, "toe_g")
        if e:
            obs.append((m, e, tx_from_meas(m, t_rx)))
    if len(obs) < 4:
        return None
    x = np.zeros(3)
    b = ephem.C * (t_rx - float(np.median([float(o[2]) for o in obs])) - 0.07)
    for it in range(iters):
        H, r = [], []
        for m, e, tx in obs:
            pos, vel, clk, _ = ephem.gps_sat(e, tx)
            pos, vel, clk, _ = ephem.gps_sat(e, tx - clk)
            th = ephem.OMEGA_E * (t_rx - tx + clk)
            sp = np.array([pos[0] * math.cos(th) + pos[1] * math.sin(th), -pos[0] * math.sin(th) + pos[1] * math.cos(th), pos[2]])
            d = sp - x
            rng = np.linalg.norm(d)
            H.append(np.r_[-d / rng, 1.0])
            r.append(float(LD(ephem.C) * (LD(t_rx) - tx)) + ephem.C * clk - rng - b)
        H, r = np.array(H), np.array(r)
        dx = np.linalg.lstsq(H, r, rcond=None)[0]
        x += dx[:3]
        b += dx[3]
        if np.linalg.norm(dx[:3]) < 1e-2:
            break
    return dict(pos=x, bias_s=b / ephem.C, resid=r - H @ dx, used=[o[0]["sv"] for o in obs])


def solve_static_batch(obs, nav, iters=15, x0=None):
    """Cold start without 4 simultaneous satellites: static receiver, continuous chip clock.
    obs: list of (epoch_ms, meas) - GPS locked measurements from many epochs (any number per
    epoch). Unknowns: position (3) + clock offset b0 (m) + drift (m/s) of the epoch time model
    t_k = t0 + (ms_k - ms0)/1000 (t0 = first measurement's transmit time + 75 ms).
    Starts from the Earth's centre unless x0 is given. Returns dict(pos, b0, drift, resid, n) or None."""
    rows = []
    for ms, m in obs:
        e = ephem.pick(nav["G"].get(m["sv"], []), m["t_hint"], "toe_g")
        if e:
            rows.append((ms, m, e, tx_from_meas(m, m["t_hint"])))
    svs = {r[1]["sv"] for r in rows}
    if len(rows) < 8 or len(svs) < 3:
        return None
    ms0 = rows[0][0]
    t0 = rows[0][3] + 0.075
    x = np.zeros(3) if x0 is None else np.array(x0, float)
    b, d = 0.0, 0.0
    for it in range(iters):
        H, r = [], []
        for ms, m, e, tx in rows:
            dt = (ms - ms0) / 1000.0
            trx = t0 + dt
            pos, vel, clk, _ = ephem.gps_sat(e, tx)
            pos, vel, clk, _ = ephem.gps_sat(e, tx - clk)
            th = ephem.OMEGA_E * (trx - tx + clk)
            sp = np.array([pos[0] * math.cos(th) + pos[1] * math.sin(th), -pos[0] * math.sin(th) + pos[1] * math.cos(th), pos[2]])
            dv = sp - x
            rng = np.linalg.norm(dv)
            H.append(np.r_[-dv / rng, 1.0, dt])
            r.append(ephem.C * (trx - tx) + ephem.C * clk - rng - b - d * dt)
        H, r = np.array(H), np.array(r)
        dx = np.linalg.lstsq(H, r, rcond=None)[0]
        x += dx[:3]
        b += dx[3]
        d += dx[4]
        if np.linalg.norm(dx[:3]) < 0.1:
            break
    res = r - H @ dx
    return dict(pos=x, b0=b, drift=d, resid=res, n=len(rows), svs=sorted(svs), cond=float(np.linalg.cond(H)))


BIT = 0.020


def tx_gps_time(m, t_hint):
    """transmit time in GPS time for a GPS (sys 1) or GLONASS (sys 0x10) FE29 measurement"""
    if m.get("sys") != 0x10:
        return tx_from_meas(m, t_hint)
    tod = LD(m["z"]) + LD(m["sub"] * 65536 + m["tx_chips"]) / LD(0.511e6)   # GLONASS time of day
    u = ephem.gpst_to_unix(t_hint)
    day0 = (u + 3 * 3600) // 86400 * 86400 - 3 * 3600                     # whole seconds: exact
    base = day0 - 3 * 3600
    base += round((u - (base + float(tod))) / 86400) * 86400
    return LD(ephem.unix_to_gpst(base)) + tod


def _eph(nav, m, t):
    if m.get("sys") == 0x10:
        return ephem.pick(nav["R"].get(m.get("slot"), []), t, "tb"), ephem.glo_sat
    return ephem.pick(nav["G"].get(m["sv"], []), t, "toe_g"), ephem.gps_sat


def resolve_bits(txs, bits=None):
    """txs: transmit times of one epoch (some may be off by whole 20 ms navigation bits although
    the chip reports them locked). Ranges of visible satellites span < 19.4 ms, so relative to a
    reference each satellite has 1-2 candidate bit corrections. Returns a list of candidate
    correction vectors (seconds), first = most compact."""
    import itertools
    ref = sorted(txs)[len(txs) // 2]           # the median: a minority of wrong channels can't
                                               # shift the whole epoch (that biases the fix by ~100 m)
    bits = bits or [BIT] * len(txs)
    win = 0.0195 if all(b == BIT for b in bits) else 0.023      # GLONASS ranges are shorter
    cands = []
    for t, bl in zip(txs, bits):
        ks = [k for k in range(-120, 121) if abs(t + k * bl - ref) <= win]
        cands.append([k * bl for k in ks] or [0.0])
    combos = list(itertools.product(*cands))[:256]
    combos.sort(key=lambda c: max(t + d for t, d in zip(txs, c)) - min(t + d for t, d in zip(txs, c)))
    return combos


def solve_abs_robust(meas, nav, t_rx):
    """solve_abs with navigation-bit ambiguity resolution (enumerated, best fit on the Earth)."""
    base = [(m, tx_from_meas(m, t_rx)) for m in meas if ephem.pick(nav["G"].get(m["sv"], []), t_rx, "toe_g")]
    if len(base) < 4:
        return None
    best = None
    for corr in resolve_bits([tx for _, tx in base]):
        ms = []
        for (m, tx), dc in zip(base, corr):
            mm = dict(m)
            mm["tx_override"] = tx + dc
            ms.append(mm)
        s = solve_abs(ms, nav, t_rx)
        if not s:
            continue
        rms = float(np.sqrt(np.mean(s["resid"] ** 2)))
        if 6.2e6 < np.linalg.norm(s["pos"]) < 6.5e6 and (best is None or rms < best[0]):
            best = (rms, s)
    return best[1] if best else None


def _ellipsoid_radius(x):
    lat = math.atan2(x[2], math.hypot(x[0], x[1]))
    a, b = ephem.WGS_A, ephem.WGS_A * (1 - ephem.WGS_F)
    c, s = math.cos(lat), math.sin(lat)
    return math.sqrt(((a * a * c) ** 2 + (b * b * s) ** 2) / ((a * c) ** 2 + (b * s) ** 2))


def solve_abs_alt(meas, nav, t_rx, h0=200.0, sigma_h=500.0, iters=30, use_dop=True):
    """Cold fix with 3+ satellites: absolute transmit times (bit ambiguities enumerated), Doppler
    range rates (they constrain the position through the line-of-sight velocity; one extra
    unknown, the receiver frequency offset) and an altitude constraint |x| = R_ellipsoid + h0.
    Several starting points on the ellipsoid below the satellites; best fit wins.
    Returns dict(pos, bias_s, drift, resid, used, rms) or None."""
    base = [(m, tx_gps_time(m, t_rx)) for m in meas if _eph(nav, m, t_rx)[0]]
    if len(base) < 3:
        return None
    med = float(np.median([float(tx) for _, tx in base]))
    base = [(m, tx) for m, tx in base if abs(float(tx) - med) < 0.2]
    if len(base) < 3:
        return None
    best = None
    for corr in resolve_bits([tx for _, tx in base], [0.010 if m.get("sys") == 0x10 else BIT for m, _ in base]):
        obs = []
        for (m, tx), dc in zip(base, corr):
            e, f = _eph(nav, m, t_rx)
            pos, vel, clk, cdr = f(e, tx + dc)
            pos, vel, clk, cdr = f(e, tx + dc - clk)
            obs.append((m, tx + dc, np.array(pos), np.array(vel), clk, cdr))
        u = np.mean([o[2] / np.linalg.norm(o[2]) for o in obs], axis=0)
        u /= np.linalg.norm(u)
        e1 = np.cross(u, [0, 0, 1.0]); e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1)
        starts = [u] + [math.cos(0.5) * u + math.sin(0.5) * (math.cos(a) * e1 + math.sin(a) * e2) for a in np.arange(0, 6.28, 1.05)]
        for st in starts:
            x = st / np.linalg.norm(st) * 6371e3
            b = ephem.C * (t_rx - float(np.median([float(o[1]) for o in obs])) - 0.075)
            dr = 0.0
            ok = True
            for it in range(iters):
                H, r, w = [], [], []
                for m, tx, sp0, sv, clk, cdr in obs:
                    th = ephem.OMEGA_E * (t_rx - tx + clk)
                    R = np.array([[math.cos(th), math.sin(th), 0], [-math.sin(th), math.cos(th), 0], [0, 0, 1.0]])
                    sp, svv = R @ sp0, R @ sv
                    d = sp - x
                    rng = np.linalg.norm(d)
                    uu = d / rng
                    H.append(np.r_[-uu, 1.0, 0.0]); r.append(float(LD(ephem.C) * (LD(t_rx) - tx)) + ephem.C * clk - rng - b); w.append(1 / 30.0)
                    if use_dop:
                        lam = ephem.C / (1602e6 + 0.5625e6 * (m["sv"] - 8)) if m.get("sys") == 0x10 else LAM
                        prr = -m["dop"] * lam + ephem.C * cdr
                        rr = uu @ svv
                        g = -(svv - uu * rr) / rng            # d(range rate)/dx
                        H.append(np.r_[g, 0.0, 1.0]); r.append(prr - rr - dr); w.append(1 / 2.0)
                nx = np.linalg.norm(x)
                H.append(np.r_[x / nx, 0.0, 0.0]); r.append(_ellipsoid_radius(x) + h0 - nx); w.append(1 / sigma_h)
                W = np.diag(w)
                H, r = np.array(H), np.array(r)
                dx = np.linalg.lstsq(W @ H, W @ r, rcond=None)[0]
                step = np.linalg.norm(dx[:3])
                if step > 3e6:
                    dx[:3] *= 3e6 / step
                x += dx[:3]; b += dx[3]; dr += dx[4]
                if np.linalg.norm(x) < 1e6:
                    ok = False
                    break
                if step < 0.1:
                    break
            if not ok:
                continue
            res = r - H @ dx
            pr = res[0:-1:2] if use_dop else res[:-1]
            score = float(np.sqrt(np.mean((W @ res) ** 2)))
            if 6.3e6 < np.linalg.norm(x) < 6.45e6 and (best is None or score < best[0]):
                best = (score, dict(pos=x.copy(), bias_s=b / ephem.C, drift=dr, resid=pr, used=[o[0]["sv"] for o in obs], rms=score))
    return best[1] if best else None


def solve_batch2(epochs, nav, h0=200.0, sigma_h=300.0, iters=40, step=4, max_pr_rms=35.0):
    """robust wrapper around _solve_batch2: a corrupted (single-decoded) orbit either spoils the
    fit (drop the satellite with the largest residuals) or prevents convergence (leave-one-out);
    always keeps 3+ satellites. Correct fits have ~20-25 m pseudorange rms."""
    keys = sorted({(m.get("sys", 1), m["sv"]) for _, ml in epochs for m in ml})

    def run(excl):
        eps = [(ms, [m for m in ml if (m.get("sys", 1), m["sv"]) not in excl]) for ms, ml in epochs]
        s = _solve_batch2(eps, nav, h0, sigma_h, iters, step)
        if s:
            s["excluded"] = sorted(excl)
        return s

    s = run(set())
    if s and s["pr_rms"] <= max_pr_rms:
        return s
    cands = []
    if s and len(s["svs"]) > 3:
        cands.append({max(s["sv_rms"], key=s["sv_rms"].get)})
    if len(keys) > 3:
        cands += [{k} for k in keys]
    for excl in cands:
        s2 = run(excl)
        if s2 and s2["pr_rms"] <= max_pr_rms and len(s2["svs"]) >= 3:
            return s2
    return s if s and s["pr_rms"] <= max_pr_rms else None


def _solve_batch2(epochs, nav, h0=200.0, sigma_h=300.0, iters=40, step=4):
    """Static receiver, 2+ satellites per epoch over minutes: pseudoranges + Doppler of all epochs
    with ONE continuous chip clock (offset, drift, drift rate over the chip's own epoch counter)
    and an altitude constraint. Navigation-bit ambiguities: each satellite's transmit time is made
    continuous from epoch to epoch (+1 s), then one 20/10 ms offset per satellite is enumerated.
    epochs: list of (ms, [meas]) with meas carrying t_hint (and slot for GLONASS).
    Returns dict(pos, rms, n, svs) or None."""
    import itertools
    rows = {}
    for ms, ml in epochs:
        for m in ml:
            e, f = _eph(nav, m, m["t_hint"])
            if not e:
                continue
            key = (m.get("sys", 1), m["sv"])
            tx = tx_gps_time(m, m["t_hint"])
            rows.setdefault(key, []).append([ms, m, e, f, tx])
    rows = {k: v for k, v in rows.items() if len(v) >= 10}
    if len(rows) < 2:
        return None
    for key, lst in rows.items():                    # continuity of each satellite's transmit time
        lst.sort(key=lambda r: r[0])
        bl = 0.010 if key[0] == 0x10 else BIT
        for a, b in zip(lst, lst[1:]):
            exp = a[4] + (b[0] - a[0]) / 1000.0
            b[4] += round((exp - b[4]) / bl) * bl if abs(exp - b[4]) > bl / 2 else 0.0
    keys = sorted(rows)
    ms0 = min(r[0] for v in rows.values() for r in v)
    first = {k: rows[k][0] for k in keys}
    # candidate per-satellite offsets relative to the first key (all tx must fit one 23 ms window)
    ref = keys[0]
    cand = []
    for k in keys:
        rk = first[k]
        bl = 0.010 if k[0] == 0x10 else BIT
        shift = (rk[0] - first[ref][0]) / 1000.0
        base = first[ref][4] + shift
        cand.append([n * bl for n in range(-150, 151) if abs(rk[4] + n * bl - base) <= 0.023] or [0.0])
    best = None
    minima = []
    combos = sorted(itertools.product(*cand), key=lambda c: sum(abs(x) for x in c))[:8]
    for combo in combos:
        off = dict(zip(keys, combo))
        DT, TX, SP, SV, CLK, CDR, PRR, SVK = [], [], [], [], [], [], [], []
        for ki, k in enumerate(keys):
            for ms, m, e, f, tx in rows[k][::step]:
                SVK.append(ki)
                tt = tx + off[k]
                pos, vel, clk, cdr = f(e, tt)
                pos, vel, clk, cdr = f(e, tt - clk)
                th = ephem.OMEGA_E * 0.072                       # Sagnac, nominal travel time
                R = np.array([[math.cos(th), math.sin(th), 0], [-math.sin(th), math.cos(th), 0], [0, 0, 1.0]])
                lam = ephem.C / (1602e6 + 0.5625e6 * (m["sv"] - 8)) if k[0] == 0x10 else LAM
                DT.append((ms - ms0) / 1000.0); TX.append(tt); SP.append(R @ np.array(pos)); SV.append(R @ np.array(vel))
                CLK.append(clk); CDR.append(cdr); PRR.append(-m["dop"] * lam)
        DT, SP, SV, CLK, CDR, PRR = map(np.array, (DT, SP, SV, CLK, CDR, PRR))
        TX = np.array(TX, dtype=LD)
        n = len(DT)
        t0 = TX.min() + LD(0.072)
        PRm = (LD(ephem.C) * (t0 + DT.astype(LD) - TX)).astype(float) + ephem.C * CLK   # measured pseudoranges
        PRRm = PRR + ephem.C * CDR
        u = (SP / np.linalg.norm(SP, axis=1)[:, None]).mean(axis=0); u /= np.linalg.norm(u)
        e1 = np.cross(u, [0, 0, 1.0]); e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1)
        for st in [u] + [math.cos(0.45) * u + math.sin(0.45) * (math.cos(a) * e1 + math.sin(a) * e2) for a in np.arange(0, 6.28, 1.05)]:
            x = st / np.linalg.norm(st) * 6371e3
            q = np.zeros(4)
            ok = True
            for it in range(iters):
                D = SP - x
                rng = np.linalg.norm(D, axis=1)
                U = D / rng[:, None]
                rr = np.einsum("ij,ij->i", U, SV)
                G = -(SV - U * rr[:, None]) / rng[:, None]
                Hp = np.c_[-U, np.ones(n), DT, DT * DT / 2, np.zeros(n)]
                rp = PRm - rng - (q[0] + q[1] * DT + q[2] * DT * DT / 2)
                Hd = np.c_[G, np.zeros(n), np.ones(n), DT, np.ones(n)]
                rd = PRRm - rr - (q[1] + q[2] * DT + q[3])
                nx = np.linalg.norm(x)
                Ha = np.r_[x / nx, 0, 0, 0, 0][None, :]
                ra = np.array([_ellipsoid_radius(x) + h0 - nx])
                H = np.vstack([Hp / 30.0, Hd / 1.0, Ha / sigma_h])
                r = np.r_[rp / 30.0, rd / 1.0, ra / sigma_h]
                dx = np.linalg.lstsq(H, r, rcond=None)[0]
                stp = np.linalg.norm(dx[:3])
                if stp > 2e6:
                    dx[:3] *= 2e6 / stp
                x += dx[:3]; q += dx[3:]
                if np.linalg.norm(x) < 1e6:
                    ok = False
                    break
                if stp < 0.5:
                    break
            if not ok or not 6.3e6 < np.linalg.norm(x) < 6.45e6:
                continue
            res = r - H @ dx
            score = float(np.sqrt(np.mean(res ** 2)))
            minima.append((score, x.copy()))
            if best is None or score < best[0]:
                prres = res[:n] * 30.0
                sv_rms = {k: float(np.sqrt(np.mean(prres[np.array(SVK) == i] ** 2))) for i, k in enumerate(keys) if np.any(np.array(SVK) == i)}
                best = (score, dict(pos=x.copy(), rms=score, n=n, svs=keys, pr_rms=float(np.sqrt(np.mean(prres ** 2))), sv_rms=sv_rms))
    if best:
        far = [sc for sc, x in minima if np.linalg.norm(x - best[1]["pos"]) > 100e3]
        best[1]["second"] = min(far) if far else float("inf")      # best distinct alternative
        best[1]["span"] = (max(r[0] for v in rows.values() for r in v) - ms0) / 1000.0
    return best[1] if best else None
