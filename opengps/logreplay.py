#!/usr/bin/python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Re-run the engine offline on a recorded session: the chip messages of an `ogps.py --log` file
are fed to ogps.Daemon in order, with the clock set to their timestamps (commands to the chip are
discarded). Use it to test solver changes on real measurements, or to debug a submitted log.

usage: logreplay.py LOG [--nav BRDC.rnx.gz] [--pos lat,lon,h | --cold] [--ref lat,lon,h] [-v]

The chip's behaviour in the log came from the commands of the session that recorded it, so this
reproduces the solvers and the bookkeeping, not the tracking.
"""
import argparse, contextlib, io, os, sys, tempfile, time, types

NOW = [0.0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("--nav", help="RINEX broadcast nav file (default: orbits only from the satellites)")
    ap.add_argument("--pos", help="start position lat,lon,h (warm); default: cold start")
    ap.add_argument("--cold", action="store_true")
    ap.add_argument("--ref", help="reference position lat,lon,h: print the error statistics")
    ap.add_argument("--fe10", default="0,0.475")
    ap.add_argument("-v", action="store_true", help="print all engine output")
    a = ap.parse_args()

    time.time = lambda: NOW[0]                  # the engine's clock follows the log
    time.monotonic = lambda: NOW[0]
    import numpy as np
    import ephem, navmsg, ogps, stpgps

    frames = []
    for line in open(a.log):
        x = line.split()
        if len(x) == 4 and x[1] == "c" and len(x[2]) == 4:
            try:
                frames.append((float(x[0]), int(x[2], 16), bytes.fromhex(x[3])))
            except ValueError:
                pass
    if not frames:
        sys.exit("no chip messages in %s" % a.log)

    class Chip:
        i = 0

        def recv(self, timeout):
            if self.i >= len(frames):
                raise EOFError
            t, ty, p = frames[self.i]
            self.i += 1
            NOW[0] = t
            return types.SimpleNamespace(mtype=ty, ok=True,
                                         payload=p if ty in ogps.PLAIN else stpgps.scramble(p))

        def send(self, *args, **kw):
            pass

    NOW[0] = frames[0][0]
    tmp = tempfile.mkdtemp(prefix="logreplay-")
    nav = ephem.load_rinex(a.nav) if a.nav else {"G": {}, "R": {}}
    cold = a.cold or not a.pos
    lat, lon, h = (float(v) for v in a.pos.split(",")) if a.pos else (0.0, 0.0, 0.0)
    d = ogps.Daemon(Chip(), nav, ephem.geo2ecef(lat, lon, h), open(os.devnull, "w"))
    d.store = navmsg.NavStore(nav, os.path.join(tmp, "ephemeris.json"))
    d.nmea, d.state, d.glonass, d.cold = None, os.path.join(tmp, "last.json"), False, cold
    d.fe10_c0, d.fe10_c1 = (float(v) for v in a.fe10.split(","))
    if not cold:
        d.prior = np.array(ephem.geo2ecef(lat, lon, h))
    out = io.StringIO()
    t0 = time.perf_counter()
    with contextlib.redirect_stdout(out):
        try:
            d.run(0)
        except EOFError:
            pass
    for line in out.getvalue().splitlines():
        if a.v or any(k in line for k in ("COLD FIX", "TIME", "EPH", "batch: candidate", "back to")):
            print(line)
    print("%d chip messages, %d fixes (%.1f s)" % (len(frames), len(d.fixes), time.perf_counter() - t0))
    if d.fixes and a.ref:
        ref = np.array(ephem.geo2ecef(*(float(v) for v in a.ref.split(","))))
        rh = float(a.ref.split(",")[2])
        e3 = [np.linalg.norm(np.array(ephem.geo2ecef(*f)) - ref) for f in d.fixes]
        eh = [np.linalg.norm(np.array(ephem.geo2ecef(f[0], f[1], rh)) - ref) for f in d.fixes]
        print("error to --ref: 3D median %.1f m (90%%: %.1f m), horizontal median %.1f m (90%%: %.1f m)"
              % (np.median(e3), np.percentile(e3, 90), np.median(eh), np.percentile(eh, 90)))


if __name__ == "__main__":
    main()
