# SPDX-License-Identifier: GPL-2.0-or-later
"""MediaTek connsys GPS (MT6755 / MT6625L) host <-> GPS-core message layer, as seen on /dev/stpgps.

Reverse engineered from stock mnld traffic (tools/gps/opengps/README.md):

  frame   = AA F0 | len:u16le | type:u16le | payload | cksum:u16le | AA 0F
  len     = low 12 bits: bytes of type + payload + cksum; top 4 bits: rolling sequence
            number (chip -> host only; the host sends 0)
  cksum   = 16-bit sum of the unescaped bytes from len up to the end of payload
  escape  = inside a frame 0xAA is sent as DE DF and 0xDE as DE E0
  type    = 0xFExx
  payload = host -> chip: plain; chip -> host: scrambled, plain = ror3(byte ^ 0xCA)
"""

ESC = 0xDE
_UNESC = {0xDF: 0xAA, 0xE0: 0xDE}
_ESC = {0xAA: b"\xde\xdf", 0xDE: b"\xde\xe0"}

# names: guesses from traffic (direction, rate, size); refine as fields are understood
NAMES = {
    0xFE00: "h.tick", 0xFE01: "h.time", 0xFE02: "h.chan", 0xFE03: "h.cfg3", 0xFE04: "h.fe04",
    0xFE05: "h.ctl", 0xFE08: "h.config", 0xFE0A: "h.fe0a", 0xFE0C: "h.fwchunk", 0xFE10: "h.fe10",
    0xFE1C: "h.sync",
    0xFE28: "c.fe28", 0xFE29: "c.meas", 0xFE2A: "c.fe2a", 0xFE2B: "c.fe2b", 0xFE2C: "c.epoch",
    0xFE30: "c.clock", 0xFE32: "c.ack", 0xFE3A: "c.fe3a", 0xFE3B: "c.fe3b", 0xFE42: "c.ready",
}


def ror3(v):
    return ((v >> 3) | (v << 5)) & 0xFF


def rol3(v):
    return ((v << 3) | (v >> 5)) & 0xFF


def descramble(p):
    return bytes(ror3(c ^ 0xCA) for c in p)


def scramble(p):
    return bytes(rol3(c) ^ 0xCA for c in p)


def unescape(s):
    o = bytearray()
    i = 0
    while i < len(s):
        if s[i] == ESC and i + 1 < len(s) and s[i + 1] in _UNESC:
            o.append(_UNESC[s[i + 1]])
            i += 2
        else:
            o.append(s[i])
            i += 1
    return bytes(o)


def encode(mtype, payload, seq=0):
    """Build a host -> chip frame."""
    n = 2 + len(payload) + 2
    body = bytes([n & 0xFF, ((n >> 8) & 0x0F) | (seq << 4)]) + mtype.to_bytes(2, "little") + payload
    body += (sum(body) & 0xFFFF).to_bytes(2, "little")
    return b"\xaa\xf0" + b"".join(_ESC.get(c, bytes([c])) for c in body) + b"\xaa\x0f"


class Frame:
    __slots__ = ("mtype", "seq", "payload", "ok")

    def __init__(self, mtype, seq, payload, ok):
        self.mtype, self.seq, self.payload, self.ok = mtype, seq, payload, ok

    @property
    def name(self):
        return NAMES.get(self.mtype, "%04x" % self.mtype)


class Deframer:
    """Incremental frame parser for one direction of the byte stream."""

    def __init__(self):
        self.buf = bytearray()
        self.junk = 0

    def feed(self, data):
        self.buf += data
        out = []
        while True:
            i = self.buf.find(b"\xaa\xf0")
            if i < 0:
                keep = 1 if self.buf[-1:] == b"\xaa" else 0
                self.junk += len(self.buf) - keep
                del self.buf[: len(self.buf) - keep]
                break
            if i:
                self.junk += i
                del self.buf[:i]
            k = self.buf.find(b"\xaa\x0f", 2)
            if k < 0:
                break
            u = unescape(bytes(self.buf[2:k]))
            del self.buf[: k + 2]
            if len(u) < 6:
                self.junk += len(u)
                continue
            ln = (u[0] | u[1] << 8) & 0x0FFF
            ok = len(u) == 2 + ln and (sum(u[:-2]) & 0xFFFF) == (u[-2] | u[-1] << 8)
            out.append(Frame(u[2] | u[3] << 8, u[1] >> 4, u[4:-2], ok))
        return out
