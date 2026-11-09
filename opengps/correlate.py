#!/usr/bin/python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Field discovery: line up FE02 (host->chip per-SV command) and FE29 (chip->host measurement)
integer fields against predicted Doppler / range from broadcast ephemeris.

usage: correlate.py FRAMES NAV.rnx.gz DATE(YYYY-MM-DD) LAT LON H [fe02|fe29]
(LAT LON H: the receiver's approximate position, e.g. from the stock driver's fix)
FRAMES = gpsmon.py --bin output: "HH:MM:SS.us dir type seq payloadhex" (chip payloads descrambled)
"""
import sys, calendar, struct, collections, math
import ephem

frames_path, nav_path, date = sys.argv[1:4]
if len(sys.argv) < 7:
    sys.exit(__doc__)
lat, lon, h = (float(x) for x in sys.argv[4:7])
rx = ephem.geo2ecef(lat, lon, h)
nav = ephem.load_rinex(nav_path)
y, mo, d = (int(x) for x in date.split("-"))
day0 = calendar.timegm((y, mo, d, 0, 0, 0))


def tparse(ts):
    hh, mm, ss = ts.split(":")
    return day0 + int(hh) * 3600 + int(mm) * 60 + float(ss)


def glo_slot_from_k(k8, t):
    for slot, eps in nav["R"].items():
        if eps[0]["k"] + 8 == k8:
            return slot
    return None


rows = collections.defaultdict(list)
for line in open(frames_path):
    ts, dr, typ, seq, hx = (line.split() + [""])[:5]
    rows[(dr, int(typ, 16))].append((tparse(ts), bytes.fromhex(hx)))

which = sys.argv[7] if len(sys.argv) > 7 else "fe02"
key = ("h", 0xFE02) if which == "fe02" else ("c", 0xFE29)
samples = collections.defaultdict(list)       # (sys, sv) -> [(t, payload)]
for t, p in rows[key]:
    if which == "fe02":
        sv, sysb = p[0], p[1]
    else:
        sv, sysb = p[1], p[2]
    samples[(sysb, sv)].append((t, p))

for (sysb, sv), lst in sorted(samples.items()):
    if sysb == 0x01:
        sc, prn = "G", sv
    elif sysb == 0x10:
        sc, prn = "R", glo_slot_from_k(sv, lst[0][0])
    else:
        continue
    if prn is None:
        continue
    print("== %s%02d (sv byte %d) n=%d" % (sc, prn, sv, len(lst)))
    step = max(1, len(lst) // 6)
    for t, p in lst[::step][:6]:
        pr = ephem.predict(nav, sc, prn, ephem.unix_to_gpst(t), rx)
        if not pr:
            continue
        freq = ephem.L1 if sc == "G" else 1602e6 + pr["eph"]["k"] * 0.5625e6
        dop = -pr["rate"] / ephem.C * freq
        i16 = [struct.unpack_from("<h", p, o)[0] for o in range(0, len(p) - 1, 2)]
        print("  t=%s dop=%8.1f el=%4.1f | %s" % (
            __import__("time").strftime("%H:%M:%S", __import__("time").gmtime(t)), dop, pr["el"],
            " ".join("%6d" % v for v in i16[:16])))
