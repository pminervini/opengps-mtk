# SPDX-License-Identifier: GPL-2.0-or-later
"""GPS LNAV navigation message from FE2B (descrambled): subframes -> broadcast ephemeris.

FE2B layout (46 bytes): [2] sv, [4] subframe counter, [6:46] 10 words x 4 bytes. The 24 data
bits of a word are bytes 0, 3, 2 (MSB..LSB); byte 1 is something else (not data). Word 1 starts
with the 0x8B preamble. Data bits arrive already corrected for the D30* inversion. Decoding follows IS-GPS-200 (subframes 1-3 = clock + ephemeris).
"""
import math

PI = 3.1415926535898


def words(p):
    return [p[6 + 4 * i] << 16 | p[9 + 4 * i] << 8 | p[8 + 4 * i] for i in range(10)]


def bits(ws, start, n):
    """start = 1-based bit index into the subframe's 24-bit words (word k covers 24(k-1)+1..24k)."""
    v = 0
    for i in range(n):
        b = start - 1 + i
        v = (v << 1) | ((ws[b // 24] >> (23 - b % 24)) & 1)
    return v


def sbits(ws, start, n):
    v = bits(ws, start, n)
    return v - (1 << n) if v >> (n - 1) else v


def cat(ws, a, na, b, nb, signed=False):
    v = (bits(ws, a, na) << nb) | bits(ws, b, nb)
    n = na + nb
    return v - (1 << n) if signed and v >> (n - 1) else v


def subframe_id(ws):
    return (ws[1] >> 2) & 7


def tow(ws):
    return (ws[1] >> 7) * 6          # TOW count of the *next* subframe start, seconds


def decode(sf):
    """sf: {1: words, 2: words, 3: words} -> ephemeris dict compatible with ephem.gps_sat"""
    s1, s2, s3 = sf[1], sf[2], sf[3]
    e = {}
    e["week"] = bits(s1, 49, 10)                       # mod 1024
    e["acc"] = bits(s1, 61, 4)
    e["health"] = bits(s1, 65, 6)
    e["iodc"] = (bits(s1, 71, 2) << 8) | bits(s1, 169, 8)
    e["tgd"] = sbits(s1, 161, 8) * 2 ** -31
    e["toc"] = bits(s1, 177, 16) * 16
    e["af2"] = sbits(s1, 193, 8) * 2 ** -55
    e["af1"] = sbits(s1, 201, 16) * 2 ** -43
    e["af0"] = sbits(s1, 217, 22) * 2 ** -31
    e["iode"] = bits(s2, 49, 8)
    e["crs"] = sbits(s2, 57, 16) * 2 ** -5
    e["dn"] = sbits(s2, 73, 16) * 2 ** -43 * PI
    e["m0"] = cat(s2, 89, 8, 97, 24, True) * 2 ** -31 * PI
    e["cuc"] = sbits(s2, 121, 16) * 2 ** -29
    e["e"] = cat(s2, 137, 8, 145, 24) * 2 ** -33
    e["cus"] = sbits(s2, 169, 16) * 2 ** -29
    e["sqrta"] = cat(s2, 185, 8, 193, 24) * 2 ** -19
    e["toe"] = bits(s2, 217, 16) * 16
    e["cic"] = sbits(s3, 49, 16) * 2 ** -29
    e["omega0"] = cat(s3, 65, 8, 73, 24, True) * 2 ** -31 * PI
    e["cis"] = sbits(s3, 97, 16) * 2 ** -29
    e["i0"] = cat(s3, 113, 8, 121, 24, True) * 2 ** -31 * PI
    e["crc"] = sbits(s3, 145, 16) * 2 ** -5
    e["w"] = cat(s3, 161, 8, 169, 24, True) * 2 ** -31 * PI
    e["omegadot"] = sbits(s3, 193, 24) * 2 ** -43 * PI
    e["iode3"] = bits(s3, 217, 8)
    e["idot"] = sbits(s3, 225, 14) * 2 ** -43 * PI
    return e


class NavStore:
    """Collects FE2B subframes per satellite; complete, consistent subframes 1-3 (same IODE/IODC)
    become an ephemeris in nav["G"][sv] (format of ephem.load_rinex). Persisted as JSON."""

    def __init__(self, nav, path=None):
        self.nav, self.path, self.sf = nav, path, {}
        self.sky = self.load(path)["G"] if path else {}

    def feed(self, p, gpst_now):
        ws = words(p)
        if ws[0] >> 16 != 0x8B:
            return None
        sid = subframe_id(ws)
        if sid not in (1, 2, 3):
            return None
        sv = p[2]
        d = self.sf.setdefault(sv, {})
        d[sid] = ws
        if not all(k in d for k in (1, 2, 3)):
            return None
        e = decode(d)
        if not (e["iode"] == e["iode3"] == (e["iodc"] & 0xFF)):
            return None
        # plausibility (the chip passes words with bit errors through; no parity is available)
        if not (5100 < e["sqrta"] < 5200 and 0 <= e["e"] < 0.05 and 0.85 < e["i0"] < 1.05
                and e["health"] == 0 and abs(e["af0"]) < 1e-3):
            return None
        wk = 2048 + e["week"]                    # 10-bit week, era 2019-04-07 .. 2038-11-20
        # toe can be in the next week near a rollover: keep toe within half a week of the TOW
        tow_now = (tow(d[1]) if 1 in d else gpst_now) % 604800
        if e["toe"] - tow_now > 302400:
            wk -= 1
        elif tow_now - e["toe"] > 302400:
            wk += 1
        e.update(prn=sv, toe_g=wk * 604800 + e["toe"], toc=wk * 604800 + e["toc"])
        # no parity from the chip: accept only an ephemeris decoded twice with identical fields
        # (subframes repeat every 30 s; bit errors do not repeat the same way)
        key = tuple(round(e[k], 15) if isinstance(e[k], float) else e[k] for k in sorted(e))
        pend = getattr(self, "pending", {})
        self.pending = pend
        if pend.get(sv) != key:
            pend[sv] = key
            self.tentative = getattr(self, "tentative", {})
            self.tentative[sv] = e                  # single decode: only for redundant (batch) solving
            d.clear()
            return None
        sky = self.sky.setdefault(sv, [])
        if any(int(x["iode"]) == e["iode"] and abs(x["toe_g"] - e["toe_g"]) < 1 for x in sky):
            return None
        sky.append(e)
        sky[:] = sky[-4:]
        lst = self.nav["G"].setdefault(sv, [])
        if not any(int(x["iode"]) == e["iode"] and abs(x["toe_g"] - e["toe_g"]) < 1 for x in lst):
            lst.append(e)
        self.save()
        return e

    def save(self):
        if self.path:
            import json
            json.dump({str(k): v for k, v in self.sky.items()}, open(self.path, "w"))

    @staticmethod
    def load(path):
        import json
        try:
            g = json.load(open(path))
        except (OSError, ValueError):
            return {"G": {}, "R": {}}
        return {"G": {int(k): v for k, v in g.items()}, "R": {}}


def page18(ws):
    """Subframe 4 page 18 (SV ID 56): Klobuchar ionosphere + UTC parameters, or None."""
    if subframe_id(ws) != 4 or bits(ws, 51, 6) != 56:
        return None
    pi = PI
    return dict(
        alpha=[sbits(ws, 57, 8) * 2 ** -30, sbits(ws, 65, 8) * 2 ** -27 / pi,
               sbits(ws, 73, 8) * 2 ** -24 / pi ** 2, sbits(ws, 81, 8) * 2 ** -24 / pi ** 3],
        beta=[sbits(ws, 89, 8) * 2 ** 11, sbits(ws, 97, 8) * 2 ** 14 / pi,
              sbits(ws, 105, 8) * 2 ** 16 / pi ** 2, sbits(ws, 113, 8) * 2 ** 16 / pi ** 3],
        a1=sbits(ws, 121, 24) * 2 ** -50, a0=cat(ws, 145, 24, 169, 8, True) * 2 ** -30,
        tot=bits(ws, 177, 8) * 2 ** 12, wnt=bits(ws, 185, 8), dtls=sbits(ws, 193, 8),
        wnlsf=bits(ws, 201, 8), dn=bits(ws, 209, 8), dtlsf=sbits(ws, 217, 8))


def klobuchar(ion, lat, lon, az, el, tow):
    """IS-GPS-200 single-frequency ionospheric delay (m). lat/lon/az/el in degrees."""
    a, b = ion["alpha"], ion["beta"]
    E = el / 180.0
    psi = 0.0137 / (E + 0.11) - 0.022
    phi_i = lat / 180.0 + psi * math.cos(math.radians(az))
    phi_i = max(-0.416, min(0.416, phi_i))
    lam_i = lon / 180.0 + psi * math.sin(math.radians(az)) / math.cos(phi_i * math.pi)
    phi_m = phi_i + 0.064 * math.cos((lam_i - 1.617) * math.pi)
    t = (4.32e4 * lam_i + tow) % 86400
    F = 1.0 + 16.0 * (0.53 - E) ** 3
    amp = max(0.0, sum(a[n] * phi_m ** n for n in range(4)))
    per = max(72000.0, sum(b[n] * phi_m ** n for n in range(4)))
    x = 2 * math.pi * (t - 50400) / per
    d = F * (5e-9 + (amp * (1 - x * x / 2 + x ** 4 / 24) if abs(x) < 1.57 else 0))
    return d * 299792458.0


def glo_data80(p):
    """FE3B (descrambled, 14 B): [0] frequency channel + 8, [2]=[3] string number, [4:14] five
    little-endian u16 words = GLONASS string bits 80..1 (MSB first; data 80..9 + Hamming 8..1)."""
    v = 0
    for i in range(4, 14, 2):
        v = (v << 16) | (p[i] | p[i + 1] << 8)
    return v


def _gb(v, hi, lo):
    return (v >> (lo - 1)) & ((1 << (hi - lo + 1)) - 1)


def _gsm(v, hi, lo):
    n = hi - lo + 1
    x = _gb(v, hi, lo)
    return -(x & ((1 << (n - 1)) - 1)) if x >> (n - 1) else x


def glo_decode(s):
    """strings {1..4: data80} -> GLONASS ephemeris (km, km/s, km/s^2 -> m), ICD L1/L2 5.1"""
    s1, s2, s3, s4 = s[1], s[2], s[3], s[4]
    return dict(
        tk=(_gb(s1, 76, 72), _gb(s1, 71, 66), _gb(s1, 65, 65) * 30),
        x=_gsm(s1, 35, 9) * 2 ** -11 * 1e3, vx=_gsm(s1, 64, 41) * 2 ** -20 * 1e3, ax=_gsm(s1, 40, 36) * 2 ** -30 * 1e3,
        health=_gb(s2, 80, 78) >> 2, tb_min=_gb(s2, 76, 70) * 15,
        y=_gsm(s2, 35, 9) * 2 ** -11 * 1e3, vy=_gsm(s2, 64, 41) * 2 ** -20 * 1e3, ay=_gsm(s2, 40, 36) * 2 ** -30 * 1e3,
        gamman=_gsm(s3, 79, 69) * 2 ** -40,
        z=_gsm(s3, 35, 9) * 2 ** -11 * 1e3, vz=_gsm(s3, 64, 41) * 2 ** -20 * 1e3, az=_gsm(s3, 40, 36) * 2 ** -30 * 1e3,
        taun=_gsm(s4, 80, 59) * 2 ** -30, en=_gb(s4, 53, 49), nt=_gb(s4, 26, 16), slot=_gb(s4, 15, 11))


# GLONASS string Hamming check masks C1..C7, taken from RTKLIB (rcvraw.c, test_glostr()),
# Copyright (c) 2007-2020 T. Takasu, BSD 2-clause license. The 8th row is C_sigma as in the
# GLONASS ICD (parity of all 85 bits, check bits b8..b1 included); RTKLIB's row covers only bits
# 85..9, which on real strings failed about half of the error-free ones.
_GLO_HAM = [[0x55, 0x55, 0x5A, 0xAA, 0xAA, 0xAA, 0xB5, 0x55, 0x6A, 0xD8, 0x08],
            [0x66, 0x66, 0x6C, 0xCC, 0xCC, 0xCC, 0xD9, 0x99, 0xB3, 0x68, 0x10],
            [0x87, 0x87, 0x8F, 0x0F, 0x0F, 0x0F, 0x1E, 0x1E, 0x3C, 0x70, 0x20],
            [0x07, 0xF8, 0x0F, 0xF0, 0x0F, 0xF0, 0x1F, 0xE0, 0x3F, 0x80, 0x40],
            [0xF8, 0x00, 0x0F, 0xFF, 0xF0, 0x00, 0x1F, 0xFF, 0xC0, 0x00, 0x80],
            [0x00, 0x00, 0x0F, 0xFF, 0xFF, 0xFF, 0xE0, 0x00, 0x00, 0x01, 0x00],
            [0xFF, 0xFF, 0xF0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x00],
            [0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xF8]]


def glo_hamming_ok(m, data80):
    """GLONASS string check (ICD extended Hamming code, all eight checks zero): 85-bit string =
    idle 0, string number (4 bits), bits 80..1; left aligned in 11 bytes. Measured on 447 real
    strings: every string matching the repeated (majority) content passes."""
    v = ((m & 15) << 80 | data80) << 3          # 85 bits -> 88, idle bit is the MSB (0)
    buf = v.to_bytes(11, "big")
    return all(sum(bin(buf[j] & mk[j]).count("1") for j in range(11)) % 2 == 0 for mk in _GLO_HAM)


class GloStore:
    """GLONASS strings (FE3B) -> ephemeris in nav["R"][slot] (format of ephem.load_rinex), accepted
    only when decoded twice identically. Also maps the frequency channel k to the slot."""

    def __init__(self, nav):
        self.nav, self.str, self.pending, self.k2slot = nav, {}, {}, {}

    def feed(self, p, gpst_now):
        ch, m = p[0], p[2]
        if not 1 <= m <= 15 or not glo_hamming_ok(m, glo_data80(p)):   # bit errors: drop the string
            return None
        d = self.str.setdefault(ch, {})
        if m == 1:
            d.clear()
        d[m] = glo_data80(p)
        if not all(k in d for k in (1, 2, 3, 4)):
            return None
        e = glo_decode(d)
        d.clear()
        r = (e["x"] ** 2 + e["y"] ** 2 + e["z"] ** 2) ** 0.5
        if not (25.0e6 < r < 26.2e6 and 1 <= e["slot"] <= 24 and e["health"] == 0):
            return None
        # identical twice = accepted; tk (frame time in string 1) changes every 30 s, so not compared
        key = tuple(sorted((k, round(v, 9) if isinstance(v, float) else v) for k, v in e.items() if k != "tk"))
        if self.pending.get(ch) != key:
            self.pending[ch] = key
            return None
        # tb (GLONASS time of day, minutes) -> GPS seconds, nearest to now
        utc_now = gpst_now - 18 + 315964800                     # GPS epoch in unix seconds
        day0 = (utc_now + 3 * 3600) // 86400 * 86400 - 3 * 3600   # start of the GLONASS (Moscow) day, unix
        tb_unix = day0 + e["tb_min"] * 60                          # (the UTC+3 offset is in day0 already)
        if tb_unix - utc_now > 43200:
            tb_unix -= 86400
        elif utc_now - tb_unix > 43200:
            tb_unix += 86400
        eph = dict(slot=e["slot"], k=ch - 8, tb=tb_unix - 315964800 + 18, taun=e["taun"], gamman=e["gamman"],
                   x=e["x"], y=e["y"], z=e["z"], vx=e["vx"], vy=e["vy"], vz=e["vz"],
                   ax=e["ax"], ay=e["ay"], az=e["az"], health=0)
        lst = self.nav["R"].setdefault(e["slot"], [])
        if any(abs(x["tb"] - eph["tb"]) < 1 for x in lst):
            return None
        lst.append(eph)
        self.k2slot[ch - 8] = e["slot"]
        return eph
