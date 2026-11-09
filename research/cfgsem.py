#!/usr/bin/python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Semantics of the required FE08 fields: which values start / track. Needs NAV=<RINEX nav> POS=lat,lon,h.
import sys, os, time, subprocess, re
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'opengps'))
OGPS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'opengps', 'ogps.py')
import gpscfg
from fwload import Chip, CHUNK
fw = open('/lib/firmware/mediatek/mt6755_gps.bin', 'rb').read()
def cfg_with(**kv):
    import struct
    b = bytearray(gpscfg.build())
    for k, v in kv.items():
        struct.pack_into('<H', b, 2 * int(k[1:]), v)
    return bytes(b)
def starts(cfg):
    c = Chip()
    try:
        c.send(0xFE08, cfg); c.wait(0xFE31)
        for i in range((len(fw)+CHUNK-1)//CHUNK):
            idx=(i+1).to_bytes(2,'little'); c.send(0xFE0C, idx+fw[i*CHUNK:(i+1)*CHUNK]); c.wait(0xFE32,1.0,idx)
        c.send(0xFE05, bytes.fromhex('3200000100')); c.wait(0xFE42, 2.0); return True
    except TimeoutError: return False
    finally: os.close(c.fd); time.sleep(0.5)
def tracks(cfg):
    open('/tmp/cfgsem.hex','w').write(cfg.hex(' ')+'\n')
    out = subprocess.run(['python3',OGPS,'--nav',os.environ['NAV'],'--secs','30',
        '--config','/tmp/cfgsem.hex','--ephcache','/tmp/cs-eph.json','--pos',os.environ['POS'],'--log','/tmp/cs.log'],capture_output=True,text=True,timeout=90).stdout
    return max([int(x) for x in re.findall(r'locked=(\d+)', out)] or [0])
for k, vals in (('u8', [0x8000, 0x0004, 0x8004, 0x4000, 0x0001]), ('u9', [0x8004, 0x8000, 0x0004, 0x8005, 0x8006]),
                ('u10', [123, 1, 50, 200, 255]), ('u11', [44, 1, 22, 88, 255]), ('u30', [1, 2, 3, 0x100])):
    print(k, [(hex(v), starts(cfg_with(**{k: v}))) for v in vals], flush=True)
for v in (0, 1, 2, 3, 4, 5):
    print('u29 =', v, '-> max locked', tracks(cfg_with(u29=v)), flush=True)
