# SPDX-License-Identifier: GPL-2.0-or-later
"""Encode a synthetic ephemeris into GPS LNAV subframes 1-3 and GLONASS strings 1-4 in the chip's
FE2B / FE3B layouts, then check the decoders return the same values."""
import math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim  # noqa: E402
import navmsg  # noqa: E402

PI = navmsg.PI


def put(ws, start, n, val):
    val &= (1 << n) - 1
    for i in range(n):
        b = start - 1 + i
        bit = (val >> (n - 1 - i)) & 1
        w, s = b // 24, 23 - b % 24
        ws[w] = (ws[w] & ~(1 << s)) | (bit << s)


def subframes(e, tow_next):
    s1, s2, s3 = [0] * 10, [0] * 10, [0] * 10
    for sid, s in ((1, s1), (2, s2), (3, s3)):
        put(s, 1, 8, 0x8B)
        put(s, 25, 17, (tow_next + 6 * (sid - 1)) // 6)
        put(s, 44, 3, sid)
    put(s1, 49, 10, e["week"] % 1024); put(s1, 65, 6, 0)
    put(s1, 71, 2, e["iodc"] >> 8); put(s1, 169, 8, e["iodc"])
    put(s1, 161, 8, round(e["tgd"] / 2 ** -31)); put(s1, 177, 16, (e["toc"] % 604800) // 16)
    put(s1, 193, 8, round(e["af2"] / 2 ** -55)); put(s1, 201, 16, round(e["af1"] / 2 ** -43))
    put(s1, 217, 22, round(e["af0"] / 2 ** -31))
    put(s2, 49, 8, e["iode"]); put(s2, 57, 16, round(e["crs"] / 2 ** -5)); put(s2, 73, 16, round(e["dn"] / PI / 2 ** -43))
    m0 = round(e["m0"] / PI / 2 ** -31); put(s2, 89, 8, m0 >> 24); put(s2, 97, 24, m0)
    put(s2, 121, 16, round(e["cuc"] / 2 ** -29))
    ee = round(e["e"] / 2 ** -33); put(s2, 137, 8, ee >> 24); put(s2, 145, 24, ee)
    put(s2, 169, 16, round(e["cus"] / 2 ** -29))
    a = round(e["sqrta"] / 2 ** -19); put(s2, 185, 8, a >> 24); put(s2, 193, 24, a)
    put(s2, 217, 16, e["toe"] // 16)
    put(s3, 49, 16, round(e["cic"] / 2 ** -29))
    o = round(e["omega0"] / PI / 2 ** -31); put(s3, 65, 8, o >> 24); put(s3, 73, 24, o)
    put(s3, 97, 16, round(e["cis"] / 2 ** -29))
    i0 = round(e["i0"] / PI / 2 ** -31); put(s3, 113, 8, i0 >> 24); put(s3, 121, 24, i0)
    put(s3, 145, 16, round(e["crc"] / 2 ** -5))
    w = round(e["w"] / PI / 2 ** -31); put(s3, 161, 8, w >> 24); put(s3, 169, 24, w)
    put(s3, 193, 24, round(e["omegadot"] / PI / 2 ** -43)); put(s3, 217, 8, e["iode"])
    put(s3, 225, 14, round(e["idot"] / PI / 2 ** -43))
    return s1, s2, s3


def fe2b(sv, ws):
    """chip layout: data bits of each word in bytes 0, 3, 2; byte 1 = status"""
    p = bytearray(46)
    p[2] = sv
    for i, w in enumerate(ws):
        p[6 + 4 * i], p[9 + 4 * i], p[8 + 4 * i], p[7 + 4 * i] = w >> 16, (w >> 8) & 0xFF, w & 0xFF, 0x03
    return bytes(p)


def test_gps_subframes_decode():
    e = sim.constellation()["G"][7][0]
    tow_next = int((sim.T0 % 604800) // 6 * 6)
    store = navmsg.NavStore({"G": {}, "R": {}})
    got = None
    for _ in range(2):                                     # accepted on the second identical decode
        for ws in subframes(e, tow_next):
            got = store.feed(fe2b(7, ws), sim.T0) or got
    assert got is not None
    lsb = dict(sqrta=2 ** -19, e=2 ** -33, m0=PI * 2 ** -31, omega0=PI * 2 ** -31, i0=PI * 2 ** -31,
               w=PI * 2 ** -31, dn=PI * 2 ** -43, omegadot=PI * 2 ** -43, idot=PI * 2 ** -43, crs=2 ** -5,
               crc=2 ** -5, cuc=2 ** -29, cus=2 ** -29, cic=2 ** -29, cis=2 ** -29, af0=2 ** -31,
               af1=2 ** -43, tgd=2 ** -31, toe=0, iode=0)
    for k, q in lsb.items():                               # equal to within half a field LSB
        assert abs(got[k] - e[k]) <= q / 2 + 1e-15, k
    assert got["toe_g"] == e["toe_g"]


def test_single_decode_is_only_tentative():
    e = sim.constellation()["G"][3][0]
    store = navmsg.NavStore({"G": {}, "R": {}})
    for ws in subframes(e, int((sim.T0 % 604800) // 6 * 6)):
        assert store.feed(fe2b(3, ws), sim.T0) is None
    assert 3 in store.tentative and store.tentative[3]["toc"] == e["toc"]


def sm(v, n):
    return ((1 << (n - 1)) | (-v)) if v < 0 else v


def with_hamming(m, data):
    """set check bits b8..b1 of a string (data = bits 80..9 in place) so that all checks pass"""
    data &= ~0xFF
    for cb in range(256):
        if navmsg.glo_hamming_ok(m, data | cb):
            return data | cb
    raise AssertionError("no check bits")


def test_glonass_hamming():
    import random
    rnd = random.Random(1)
    for _ in range(20):
        m = rnd.randrange(1, 16)
        d = with_hamming(m, rnd.getrandbits(80))
        assert navmsg.glo_hamming_ok(m, d)
        for bit in rnd.sample(range(80), 5):                # any single bit error is detected
            assert not navmsg.glo_hamming_ok(m, d ^ (1 << bit))
        a, b = rnd.sample(range(80), 2)                     # and any double error
        assert not navmsg.glo_hamming_ok(m, d ^ (1 << a) ^ (1 << b))


def test_glonass_strings_decode():
    x, y, z = 12345.678, -20000.25, 8000.5                # km
    vx, vy, vz = 1.234, -2.5, 3.0625                      # km/s
    s = {m: 0 for m in (1, 2, 3, 4)}

    def setb(m, hi, lo, v):
        s[m] |= (v & ((1 << (hi - lo + 1)) - 1)) << (lo - 1)
    setb(1, 76, 72, 15); setb(1, 71, 66, 46); setb(1, 65, 65, 1)
    setb(1, 35, 9, sm(round(x / 2 ** -11), 27)); setb(1, 64, 41, sm(round(vx / 2 ** -20), 24))
    setb(2, 76, 70, 63); setb(2, 35, 9, sm(round(y / 2 ** -11), 27)); setb(2, 64, 41, sm(round(vy / 2 ** -20), 24))
    setb(3, 35, 9, sm(round(z / 2 ** -11), 27)); setb(3, 64, 41, sm(round(vz / 2 ** -20), 24))
    setb(4, 80, 59, sm(round(3.6e-4 / 2 ** -30), 22)); setb(4, 15, 11, 21); setb(4, 26, 16, 1044)
    frames = {}
    s = {m: with_hamming(m, v) for m, v in s.items()}
    for m, v in s.items():                               # FE3B: 5 little-endian u16 words, MSB first
        p = bytearray(14)
        p[0], p[2], p[3] = 12, m, m
        for i in range(5):
            w = (v >> (16 * (4 - i))) & 0xFFFF
            p[4 + 2 * i], p[5 + 2 * i] = w & 0xFF, w >> 8
        assert navmsg.glo_data80(bytes(p)) == v
        frames[m] = navmsg.glo_data80(bytes(p))
    e = navmsg.glo_decode(frames)
    assert e["tk"] == (15, 46, 30) and e["tb_min"] == 63 * 15 and e["slot"] == 21 and e["nt"] == 1044
    assert math.isclose(e["x"], x * 1e3, abs_tol=1) and math.isclose(e["y"], y * 1e3, abs_tol=1)
    assert math.isclose(e["vz"], vz * 1e3, abs_tol=1e-3) and math.isclose(e["taun"], 3.6e-4, rel_tol=1e-6)
