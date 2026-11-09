#!/usr/bin/python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Open-source GPS-core bring-up for the MT6755 connsys GPS (replaces mnld's firmware download).

usage: fwload.py [--fw /lib/firmware/mediatek/mt6755_gps.bin] [--listen SECS] [--log FILE]

Sequence (from stock mnld, see README.md):
  open /dev/stpgps (powers the GPS function on via WMT), ioctls 11/10/9/7/8 (informational)
  -> FE08 config (106 bytes)              <- FE31
  -> FE0C chunk i = u16 index (1..44) + 1400 firmware bytes   <- FE32 u16 index
  -> FE05 32 00 00 01 00 (start)          <- FE42 (core alive; then 1 Hz FE42/FE2C)
With --listen the chip's messages are decoded and printed for SECS seconds. No channel
assignment is sent, so the core runs idle (that part is still mnld's job in phase 1).
The device must not be held by mnld at the same time.
"""
import argparse, fcntl, os, select, sys, time
from stpgps import Deframer, encode, descramble, NAMES
import gpscfg

HERE = os.path.dirname(os.path.abspath(__file__))
CHUNK = 1400
IOC_CO_CLOCK, IOC_RTC, IOC_EFUSE, IOC_HWVER, IOC_FWVER = 11, 10, 9, 7, 8
PLAIN = {0xFE31, 0xFE32}


class Chip:
    def __init__(self, path="/dev/stpgps"):
        self.fd = os.open(path, os.O_RDWR)
        self.rx = Deframer()
        self.q = []

    def send(self, mtype, payload=b""):
        data = encode(mtype, payload)
        for _ in range(200):                   # STP tx buffer full: wait for it to drain
            try:
                os.write(self.fd, data)
                return
            except OSError as e:
                if e.errno != 28:
                    raise
                time.sleep(0.005)
        raise OSError(28, "stpgps tx buffer stays full")

    def recv(self, timeout):
        end = time.monotonic() + timeout
        while not self.q:
            left = end - time.monotonic()
            if left <= 0 or not select.select([self.fd], [], [], left)[0]:
                return None
            self.q += self.rx.feed(os.read(self.fd, 4096))
        return self.q.pop(0)

    def wait(self, mtype, timeout=1.0, payload=None):
        end = time.monotonic() + timeout
        while True:
            f = self.recv(max(0, end - time.monotonic()))
            if f is None:
                raise TimeoutError("no %s" % NAMES.get(mtype, hex(mtype)))
            if f.mtype == mtype and (payload is None or f.payload[:len(payload)] == payload):
                return f

    def ioctl(self, cmd, arg=0):
        try:
            return fcntl.ioctl(self.fd, cmd, arg)
        except OSError as e:
            return "err %d" % e.errno


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fw", default="/lib/firmware/mediatek/mt6755_gps.bin")
    ap.add_argument("--config", help="FE08 config as hex file (default: built by gpscfg.py)")
    ap.add_argument("--minimal", action="store_true", help="only the FE08 fields the core needs")
    ap.add_argument("--listen", type=float, default=0)
    ap.add_argument("--log")
    a = ap.parse_args()
    fw = open(a.fw, "rb").read()
    cfg = bytes.fromhex("".join(l.split("#")[0] for l in open(a.config))) if a.config else gpscfg.build(a.minimal)
    t0 = time.monotonic()
    c = Chip()
    print("ioctl co_clock=%s rtc=%s" % (c.ioctl(IOC_CO_CLOCK), c.ioctl(IOC_RTC)))
    c.send(0xFE08, cfg)
    print("config  -> ack %s" % c.wait(0xFE31).payload.hex(" "))
    n = (len(fw) + CHUNK - 1) // CHUNK
    for i in range(n):
        idx = (i + 1).to_bytes(2, "little")
        c.send(0xFE0C, idx + fw[i * CHUNK:(i + 1) * CHUNK])
        c.wait(0xFE32, 1.0, idx)
    print("firmware %d bytes in %d chunks, all acked (%.0f ms)" % (len(fw), n, (time.monotonic() - t0) * 1e3))
    c.send(0xFE05, bytes.fromhex("3200000100"))
    f = c.wait(0xFE42, 2.0)
    print("start   -> FE42 after %.0f ms: GPS core firmware running" % ((time.monotonic() - t0) * 1e3))
    log = open(a.log, "w") if a.log else None
    end = time.monotonic() + a.listen
    while time.monotonic() < end:
        f = c.recv(end - time.monotonic())
        if not f:
            continue
        p = f.payload if f.mtype in PLAIN else descramble(f.payload)
        line = "%.3f <- %-9s seq%d len%3d %s" % (time.monotonic() - t0, f.name, f.seq, len(p), p.hex(" "))
        print(line[:200], flush=True)
        if log:
            log.write(line + "\n")
    os.close(c.fd)


if __name__ == "__main__":
    main()
