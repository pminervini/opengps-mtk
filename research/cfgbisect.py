#!/usr/bin/python3
# SPDX-License-Identifier: GPL-2.0-or-later
# Which FE08 config fields must be non-zero for the GPS core to start (FE42)? (run on the device)
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'opengps'))
from fwload import Chip, CHUNK
fw = open('/lib/firmware/mediatek/mt6755_gps.bin', 'rb').read()
import gpscfg
cfg = gpscfg.build()
def starts(c_bytes):
    c = Chip()
    try:
        c.send(0xFE08, c_bytes); c.wait(0xFE31)
        for i in range((len(fw)+CHUNK-1)//CHUNK):
            idx=(i+1).to_bytes(2,'little'); c.send(0xFE0C, idx+fw[i*CHUNK:(i+1)*CHUNK]); c.wait(0xFE32,1.0,idx)
        c.send(0xFE05, bytes.fromhex('3200000100'))
        c.wait(0xFE42, 2.0); return True
    except TimeoutError:
        return False
    finally:
        os.close(c.fd); time.sleep(0.5)
mode = sys.argv[1]
if mode == 'single':
    need = []
    for k in range(1, len(cfg)//2):
        b = bytearray(cfg); b[2*k:2*k+2] = b'\0\0'
        if bytes(b) == cfg: continue
        ok = starts(bytes(b))
        print('zero u16[%2d] (bytes %3d-%3d, was 0x%04x): %s' % (k, 2*k, 2*k+1, cfg[2*k] | cfg[2*k+1] << 8, 'starts' if ok else 'NO START'), flush=True)
        if not ok: need.append(k)
    print('needed:', need)
elif mode == 'keep':
    keep = [int(x) for x in sys.argv[2].split(',')]
    b = bytearray(len(cfg)); b[0:2] = cfg[0:2]
    for k in keep: b[2*k:2*k+2] = cfg[2*k:2*k+2]
    print('only fields', keep, '->', 'starts' if starts(bytes(b)) else 'NO START')
