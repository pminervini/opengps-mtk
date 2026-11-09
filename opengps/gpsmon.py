#!/usr/bin/python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Decode and log the host <-> GPS-core traffic of a running stock mnld (or a saved capture).

Input is strace output of mnld (-f -tt -y -xx -s 65535 -e trace=read,write), either a file or live:

  sudo strace -f -tt -y -xx -s 65535 -e trace=read,write -e signal=none -p $(pidof mnld) 2>&1 \
      | gpsmon.py -            # live
  gpsmon.py capture.strace     # offline

Options: -t TYPE[,TYPE..] only these message types (hex, e.g. fe29), -s summary only,
         --raw  print undescrambled chip payloads, --bin FILE  save frames (see stpgps.py) as
         "t dir type seq payloadhex" lines for later analysis.
"""
import argparse, collections, re, sys, time
from stpgps import Deframer, descramble, NAMES

LINE = re.compile(r'^(?:\[pid\s+(?P<p1>\d+)\]\s+|(?P<p2>\d+)\s+)?(?P<ts>\d\d:\d\d:\d\d\.\d+)\s+(?:(?P<op>read|write)\((?P<fd>\d+)<(?P<path>[^>]*)>(?P<rest>.*)|<\.\.\. (?P<rop>read|write) resumed>(?P<rrest>.*))$')
DATA = re.compile(r', "((?:\\x[0-9a-f]{2})*)"(?:\.\.\.)?, \d+\) = (\d+)')
STPGPS = "/dev/stpgps"

# chip -> host messages sent by the boot ROM before the firmware runs are not scrambled
PLAIN_CHIP = {0xFE31, 0xFE32}


def path(s):
    return re.sub(r"\\x([0-9a-f]{2})", lambda m: chr(int(m.group(1), 16)), s)


def meas(p):
    """FE29 (descrambled): per-channel tracking report. Field meanings are provisional."""
    return "ch=%2d sv=%3d sys=%02x cn0?=%2d" % (p[0], p[1], p[2], p[3])


DECODERS = {0xFE29: meas}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("-t", "--types")
    ap.add_argument("-s", "--summary", action="store_true")
    ap.add_argument("--raw", action="store_true")
    ap.add_argument("--bin")
    a = ap.parse_args()
    only = {int(x, 16) for x in a.types.split(",")} if a.types else None
    f = sys.stdin if a.src == "-" else open(a.src, errors="replace")
    out = open(a.bin, "w") if a.bin else None
    rx, tx = Deframer(), Deframer()
    pending = {}          # pid -> op of an unfinished read/write on /dev/stpgps
    stats = collections.Counter()
    bad = collections.Counter()
    for line in f:
        m = LINE.match(line.rstrip("\n"))
        if not m:
            continue
        pid, ts = m["p1"] or m["p2"] or "-", m["ts"]
        if m["op"]:
            op, fdpath, rest = m["op"], path(m["path"]), m["rest"]
            if not fdpath.endswith(STPGPS):
                continue
            if rest.strip() == "<unfinished ...>":
                pending[pid] = op
                continue
        else:
            op, rest = m["rop"], m["rrest"]
            if pending.pop(pid, None) != op:
                continue
        d = DATA.search(rest)
        if not d:
            continue
        data = bytes(int(x, 16) for x in re.findall(r"\\x([0-9a-f]{2})", d.group(1)))
        tohost = op == "read"
        for fr in (rx if tohost else tx).feed(data):
            key = ("c<" if tohost else "h>") + fr.name
            stats[key] += 1
            if not fr.ok:
                bad[key] += 1
            p = fr.payload
            if tohost and fr.mtype not in PLAIN_CHIP and not a.raw:
                p = descramble(p)
            if out:
                out.write("%s %s %04x %d %s\n" % (ts, "c" if tohost else "h", fr.mtype, fr.seq, p.hex()))
            if a.summary or (only and fr.mtype not in only):
                continue
            dec = DECODERS.get(fr.mtype)
            extra = (" | " + dec(p)) if dec and len(p) >= 4 else ""
            print("%s %s %-10s seq%-2d len%4d%s%s %s" % (
                ts, "<-" if tohost else "->", fr.name, fr.seq, len(p), "" if fr.ok else " BADCK",
                extra, p[:48].hex(" ") + (" ..." if len(p) > 48 else "")), flush=True)
    print("# frames:", dict(sorted(stats.items())), file=sys.stderr)
    if bad:
        print("# bad checksum:", dict(bad), file=sys.stderr)
    print("# unparsed bytes: to-chip %d, from-chip %d" % (tx.junk, rx.junk), file=sys.stderr)


if __name__ == "__main__":
    main()
