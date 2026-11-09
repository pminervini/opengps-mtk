#!/usr/bin/python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Replay test (phase-2 step 1/2): drive the GPS core with stock mnld's *captured* host commands.

usage: replay.py FRAMES [--shift] [--secs N]
FRAMES = gpsmon.py --bin capture of a stock session ("ts dir type seq payloadhex").
Host commands are grouped per second (one mnld epoch) and replayed one group per chip epoch:
  default : as captured, only FE00/FE01 epoch-counter fields replaced by the live counter
  --shift : also move all time fields (FE01 Z-count, FE02 transmit time Z/hi/chips) by the time
            elapsed since the capture, i.e. "modified" commands
Prints how many satellites the chip locks (state 0x40). Fully generated commands = ogps.py.
"""
import argparse, struct, time, os
import ephem, pvt, gpscfg
from fwload import Chip
from stpgps import descramble
from ogps import bringup, INIT_FE05

ap = argparse.ArgumentParser()
ap.add_argument("frames")
ap.add_argument("--shift", action="store_true")
ap.add_argument("--secs", type=float, default=60)
ap.add_argument("--day", default=time.strftime("%Y-%m-%d", time.gmtime()))
a = ap.parse_args()

import calendar
y, mo, d = (int(x) for x in a.day.split("-"))
day0 = calendar.timegm((y, mo, d, 0, 0, 0))
groups, cur, last_t = [], [], None
for line in open(a.frames):
    ts, dr, ty, seq, hx = (line.split() + [""])[:5]
    if dr != "h" or ty in ("fe0c",) or (ty == "fe08" and len(hx) > 40):
        continue
    hh, mm, ss = ts.split(":")
    t = day0 + int(hh) * 3600 + int(mm) * 60 + float(ss)
    if ty == "fe00" and cur:
        groups.append((cur[0][0], cur))
        cur = []
    cur.append((t, int(ty, 16), bytes.fromhex(hx)))
groups = groups[1:]                     # first group = bring-up FE05 block (sent by bringup)
print("captured epochs:", len(groups))

c = Chip()
bringup(c, open("/lib/firmware/mediatek/mt6755_gps.bin", "rb").read(), gpscfg.build())
for hx in INIT_FE05:
    c.send(0xFE05, bytes.fromhex(hx))
t0 = time.monotonic()
k, locked_max, locked = 0, 0, {}
dt_cap = None
while time.monotonic() - t0 < a.secs and k < len(groups):
    f = c.recv(1.0)
    if not f:
        continue
    p = descramble(f.payload) if f.mtype not in (0xFE31, 0xFE32) else f.payload
    if f.mtype == 0xFE29 and p[6] == 0x40:
        locked[(p[1], p[2])] = time.monotonic()
    if f.mtype != 0xFE2C:
        continue
    ms = (p[0] | p[1] << 8) << 16 | (p[2] | p[3] << 8)
    tcap, cmds = groups[k]
    if dt_cap is None:
        dt_cap = time.time() - tcap                       # capture -> now
    for tc, ty, pl in cmds:
        pl = bytearray(pl)
        if ty == 0xFE00:
            struct.pack_into("<HH", pl, 0, (ms + 1000) >> 16 & 0xFFFF, (ms + 1000) & 0xFFFF)
        elif ty == 0xFE01:
            struct.pack_into("<HH", pl, 6, ms >> 16 & 0xFFFF, ms & 0xFFFF)
            if a.shift:
                z = struct.unpack_from("<H", pl, 13)[0]
                struct.pack_into("<H", pl, 13, (z + int(round(dt_cap / 1.5))) & 0xFFFF)
        elif ty == 0xFE02 and a.shift and pl[1] == 0x01:
            z, hi, _, chips = struct.unpack_from("<HBBH", pl, 4)
            tx = z * 1.5 + (hi * 65536 + chips) / 1.023e6 + dt_cap
            z2 = int(tx // 1.5)
            n = int(round((tx - z2 * 1.5) * 1.023e6))
            struct.pack_into("<HBBH", pl, 4, z2 & 0xFFFF, n >> 16, 0, n & 0xFFFF)
        c.send(ty, bytes(pl))
    k += 1
    now = time.monotonic()
    nl = sum(1 for v in locked.values() if now - v < 2)
    locked_max = max(locked_max, nl)
    if k % 10 == 0:
        print("epoch %3d (capture +%.0f s): locked %d" % (k, tcap - groups[0][0], nl), flush=True)
os.close(c.fd)
print("replay%s: max satellites locked %d" % (" (time-shifted)" if a.shift else " (as captured)", locked_max))
