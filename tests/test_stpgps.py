# SPDX-License-Identifier: GPL-2.0-or-later
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "opengps"))
import stpgps  # noqa: E402


def test_scramble_roundtrip():
    data = bytes(range(256))
    assert stpgps.descramble(stpgps.scramble(data)) == data
    assert stpgps.descramble(b"\xca") == b"\x00"          # zero bytes appear as 0xCA on the wire


def test_frame_roundtrip_with_escapes():
    payload = b"\x01\xaa\x02\xde\xaa\xde\x00\xff"
    frame = stpgps.encode(0xFE29, payload)
    assert frame[:2] == b"\xaa\xf0" and frame[-2:] == b"\xaa\x0f"
    assert b"\xaa" not in frame[2:-2]                      # escaped inside the frame
    d = stpgps.Deframer()
    out = d.feed(frame[:5]) + d.feed(frame[5:])             # split across reads
    assert len(out) == 1 and out[0].ok and out[0].mtype == 0xFE29 and out[0].payload == payload


def test_bad_checksum_flagged():
    frame = bytearray(stpgps.encode(0xFE02, b"\x05\x01\x01\x00"))
    frame[-3] ^= 0x01                                       # corrupt the checksum
    out = stpgps.Deframer().feed(bytes(frame))
    assert len(out) == 1 and not out[0].ok


def test_sequence_nibble():
    frame = stpgps.encode(0xFE42, b"\x00" * 4, seq=5)
    out = stpgps.Deframer().feed(frame)
    assert out[0].seq == 5 and out[0].ok
