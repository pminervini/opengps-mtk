#!/bin/sh
# SPDX-License-Identifier: GPL-2.0-or-later
# Effect of the FE10 crystal-model offset on the measured residual Doppler. Needs NAV, POS.
: "${NAV:?set NAV=path/to/BRDC...rnx.gz}" "${POS:?set POS=lat,lon,h}"; OUT=${OUT:-/tmp/opengps-research}; mkdir -p $OUT
cd "$(dirname "$(readlink -f "$0")")/../opengps"
run() { sudo timeout 70 python3 ogps.py --nav "$NAV" --pos "$POS" --secs 40 --ephcache /tmp/x-eph.json --state /tmp/x-last.json --log /tmp/fx.log "$@" > /tmp/fx.out 2>&1
  ml=$(grep -o 'locked=[0-9]*' /tmp/fx.out | cut -d= -f2 | sort -n | tail -1); nf=$(grep -c '^FIX' /tmp/fx.out)
  dop=$(grep clock: /tmp/fx.log | tail -15 | sed 's/.*dop_off \([-0-9.]*\) Hz.*/\1/' | sort -n | awk '{a[NR]=$1} END{print a[int(NR/2)+1]}')
  echo "$* -> max_locked=$ml fixes=$nf median dop_off=$dop Hz"; }
run --fe10 416,0.475
run --fe10 0,0.475
run --fe10 816,0.475
run --fe01drift 0
