#!/bin/sh
# SPDX-License-Identifier: GPL-2.0-or-later
# coldtest.sh [MIN] [OUTDIR]: true cold start test - no orbits, no position, no crystal model,
# nothing cached. Prints the cold fix and the number of fixes.
M=${1:-20}; O=${2:-/tmp/opengps-coldtest}; D=$(dirname "$(readlink -f "$0")")
mkdir -p "$O"; rm -f "$O/eph.json" "$O/last.json"
sudo timeout $((M*60+60)) python3 "$D/ogps.py" --cold --fe10 0,0.475 --secs $((M*60)) \
  --ephcache "$O/eph.json" --state "$O/last.json" --nmea "$O/cold.nmea" --log "$O/cold.log" > "$O/cold.out" 2>&1
grep -E "COLD FIX|EPH|TIME|back to" "$O/cold.out"
echo "first fix: $(grep -m1 '^FIX' "$O/cold.out")"; echo "fixes: $(grep -c '^FIX' "$O/cold.out")"
