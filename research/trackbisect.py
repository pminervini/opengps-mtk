#!/usr/bin/python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Which FE08 fields are needed for tracking (>= 3 locked within 30 s), by bisection. Needs NAV, POS.
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'opengps'))
import os, subprocess, re, sys
OGPS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'opengps', 'ogps.py')
import gpscfg
cfg = gpscfg.build()
BASE = [8, 9, 10, 11, 30]
cand = [k for k in range(1, len(cfg)//2) if k not in BASE and cfg[2*k:2*k+2] != b'\0\0']
def tracks(keep):
    b = bytearray(len(cfg)); b[0:2] = cfg[0:2]
    for k in BASE + keep: b[2*k:2*k+2] = cfg[2*k:2*k+2]
    open('/tmp/cfg-test.hex', 'w').write(b.hex(' ') + '\n')
    out = subprocess.run(['python3', OGPS, '--nav', os.environ['NAV'],
                          '--secs', '28', '--config', '/tmp/cfg-test.hex', '--ephcache', '/tmp/eph-test.json', '--log', '/tmp/tb.log'],
                         capture_output=True, text=True, timeout=90).stdout
    mx = max([int(x) for x in re.findall(r'locked=(\d+)', out)] or [0])
    print('  keep %s -> max locked %d' % (keep, mx), flush=True)
    return mx >= 3
def find(group, others):
    """minimal subset of group needed (together with others) for tracking"""
    if tracks(others):
        return []
    if len(group) == 1:
        return group
    a, b = group[:len(group)//2], group[len(group)//2:]
    na = find(a, others + b)          # what of a is needed given all of b
    nb = find(b, others + na)         # what of b is needed given needed part of a
    return na + nb
print('candidates', cand, flush=True)
if not tracks(cand):
    print('full set does not track?!'); sys.exit(1)
need = find(cand, [])
print('NEEDED for tracking:', need, [hex(cfg[2*k] | cfg[2*k+1] << 8) for k in need])
