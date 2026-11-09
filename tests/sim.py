# SPDX-License-Identifier: GPL-2.0-or-later
"""Synthetic GPS constellation and FE29 measurements for the tests (no recorded data).

Orbits are made-up but realistic Kepler elements (GPS-like: a = 26560 km, i = 55 deg, six planes);
the receiver sits at a public location (Royal Observatory, Greenwich).
"""
import math, os, struct, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "opengps"))
import ephem, pvt  # noqa: E402

WEEK = 2444
T0 = WEEK * 604800 + 2 * 86400 + 36000.0          # GPS time used by the tests
RX = ephem.geo2ecef(51.4769, -0.0005, 70.0)       # Greenwich


def constellation(n=24):
    nav = {"G": {}, "R": {}}
    for k in range(n):
        plane, slot = k % 6, k // 6
        e = dict(prn=k + 1, sqrta=math.sqrt(26559.7e3), e=0.005 + 0.001 * slot, i0=math.radians(55.0),
                 omega0=plane * math.pi / 3, w=0.3 * slot, m0=slot * math.pi / 2 + plane * 0.35, dn=4.5e-9,
                 omegadot=-8.0e-9, idot=1e-10, cuc=1e-6, cus=5e-6, crc=200.0, crs=20.0, cic=5e-8, cis=-5e-8,
                 af0=1e-5 * (k - 12), af1=1e-12, af2=0.0, tgd=-5e-9, health=0, iode=k, iodc=k, week=WEEK,
                 toe=int((T0 - WEEK * 604800) // 16 * 16))
        e["toe_g"] = WEEK * 604800 + e["toe"]
        e["toc"] = e["toe_g"]
        nav["G"][k + 1] = [e]
    return nav


def fe29(sv, tx_sv, dop, state=0x40, cn0=35):
    """FE29 payload (descrambled) for satellite sv with SV-time transmit time tx_sv"""
    tow = float(tx_sv - pvt.LD(WEEK * 604800))      # exact week-relative time (tx_sv is extended precision)
    z = int(tow // 1.5)
    n = (tow - z * 1.5) * 1.023e6
    hi, chips = int(n) >> 16, int(n) & 0xFFFF
    frac = int((n - int(n)) * 65536)
    p = bytearray(57)
    p[0], p[1], p[2], p[3], p[5], p[6] = 0, sv, 1, cn0, 0x14, state
    struct.pack_into("<HBB", p, 23, z & 0xFFFF, hi, 0)
    struct.pack_into("<HHh", p, 27, chips, frac, int(round(dop)))
    return bytes(p)


def tx_time(p, t_rx, tropo=False):
    """true SV-time transmit time, in extended precision (absolute float64 GPS seconds only
    resolve 0.24 us = 72 m)"""
    delay = (p["range"] + (pvt.tropo(p["el"]) if tropo else 0.0)) / ephem.C
    return pvt.LD(t_rx) - pvt.LD(delay) + pvt.LD(p["clk"])


def measurements(nav, t_rx, rx=RX, clock_bias=1.234e-4, min_el=10, sats=None):
    """FE29 measurements of all satellites above min_el at receive time t_rx (true GPS time),
    with a receiver clock bias. Returns (list of parse_fe29 dicts, chip epoch time)."""
    out = []
    for prn in nav["G"]:
        if sats and prn not in sats:
            continue
        p = ephem.predict(nav, "G", prn, t_rx, rx)
        if p["el"] < min_el:
            continue
        tx_sv = tx_time(p, t_rx, tropo=True)
        dop = -p["rate"] / ephem.C * ephem.L1
        m = pvt.parse_fe29(fe29(prn, tx_sv, dop))
        m["t_hint"] = t_rx + clock_bias
        m["sys"] = 1
        out.append(m)
    return out, t_rx + clock_bias
