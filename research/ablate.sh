#!/bin/sh
# SPDX-License-Identifier: GPL-2.0-or-later
# Leave out each host message type in turn (45 s each) and report locks/fixes. Needs NAV, POS.
: "${NAV:?set NAV=path/to/BRDC...rnx.gz}" "${POS:?set POS=lat,lon,h}"; OUT=${OUT:-/tmp/opengps-research}; mkdir -p $OUT
cd "$(dirname "$(readlink -f "$0")")/../opengps"
for o in none fe10 fe01 init05 fe08 fe04 fe1c fe0529 fe03b fe10,fe01,init05,fe08,fe04,fe1c,fe0529,fe03b; do
  out=$(sudo timeout 75 python3 ogps.py --nav "$NAV" --pos "$POS" --secs 45 --debug-omit $( [ $o = none ] && echo "" || echo $o) --ephcache /tmp/abl-eph.json --state /tmp/abl-last.json --log /tmp/abl-$o.log 2>&1)
  ml=$(echo "$out" | grep -o 'locked=[0-9]*' | cut -d= -f2 | sort -n | tail -1)
  nf=$(echo "$out" | grep -c '^FIX')
  first=$(echo "$out" | grep -m1 '^FIX' | cut -d' ' -f2)
  med=$(echo "$out" | grep fixes: )
  echo "omit=$o  max_locked=$ml  fixes=$nf first=$first  $med"
done
