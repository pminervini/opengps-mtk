# opengps-mtk

An open-source GPS engine for the GPS core in MediaTek's connectivity subsystem (connsys). It
replaces MediaTek's proprietary positioning stack (`mnld` + `libmnl.so`) on mainline Linux:
it loads the firmware, drives the GPS core every second, decodes its measurements and the
navigation data from the satellites, and computes position, velocity and time. It writes NMEA,
so gpsd and everything above it works unchanged.

No Android userspace, no chroot, no binary blobs other than the GPS-core firmware, which
`fwextract.py` cuts out of the device's own `libmnl.so` (in the spirit of b43-fwcutter).

Written in Python 3 with numpy. GPL-2.0-or-later.

## Status

Developed and tested on one device: an **MT6750 tablet** (MT6755-family SoC, MT6625L RF; the
vendor `libmnl.so` comes from its Android 8 image), running mainline Linux.

What works:
- **Warm start** (approximate position and clock known): first fix in 7-20 s; then one fix per
  second from an extended Kalman filter, 6-10 satellites at a window. Fixes are within 5-20 m
  of the stock engine's on the same device.
- **Cold start** (nothing known: no position, no orbits, no almanac, uncalibrated crystal): the
  core searches the sky; GPS time and orbits are decoded from the satellites. At a window with
  only 2-3 satellites in view at a time, the first fix came after about 24 minutes and was
  103 m from a surveyed reference. The stock engine got no fix in 15 minutes at the same spot.
  With open sky and 4+ satellites, a single-epoch fix is possible (52 m median in replays) as
  soon as the orbits are decoded.
- **Assistance from the Internet** (optional): today's broadcast ephemeris from IGS (BKG open
  data) with `--download`.
- **GPS L1 C/A**; **GLONASS** decoding and solving are implemented but less tested.

Other MediaTek SoCs with the same connsys GPS core (MT67xx) are likely similar but untested. The
firmware image is selected by its SHA-256, so other chips need their hash added to
`fwextract.py`. Captures from other devices are very welcome (see [Contributing](#contributing)).

## Requirements

- A kernel driver that provides MediaTek's `/dev/stpgps` character device (the GPS channel of
  the WMT/STP connsys transport). Downstream Android kernels have it. For mainline, there is the
  out-of-tree `conn_soc` driver from
  [mt6572-connsys-modules](https://github.com/buddyjojo/mt6572-connsys-modules), which needs to
  be adapted to your SoC.
- The GPS RF front-end powered: the right PMIC supplies (on MT6750 with the MT6353 PMIC this
  includes VTCXO24) and the LNA enable GPIO, as on the stock kernel.
- Python 3 (3.9 or later) and numpy.
- The device's vendor `libmnl.so`, to extract the firmware from.

## Quick start

```sh
# 1. Firmware, once: from the device's vendor partition (mounted at /vendor, or give the path)
sudo python3 opengps/fwextract.py /path/to/vendor/lib/libmnl.so

# 2. Check that the core boots and talks (prints decoded messages for 10 s)
sudo python3 opengps/fwload.py --listen 10

# 3. Run it. --pos for a warm start; without --pos (and with no last fix saved) it starts cold.
sudo python3 opengps/ogps.py --download --pos 51.48,-0.0,50 --nmea /dev/stdout
sudo python3 opengps/ogps.py --cold --nmea /tmp/gps.nmea
```

The last fix and the learned crystal model are saved in `/var/lib/opengps/last.json`, and
ephemerides decoded from the satellites in `/var/lib/opengps/ephemeris.json`, so the next start
is warm.

### As a service, with gpsd

`systemd/opengps.service` runs `ogps.py --download --pty /run/opengps/gps0`. `--pty` creates a
pseudo-terminal and symlinks it there, so gpsd can use it like a serial GPS:

```sh
sudo cp -r opengps /usr/local/lib/opengps
sudo cp systemd/opengps.service /etc/systemd/system/
sudo systemctl enable --now opengps
# /etc/default/gpsd:  DEVICES="/run/opengps/gps0"  GPSD_OPTIONS="-n"
sudo systemctl restart gpsd
gpspipe -w | head
```

### Useful options

| option | meaning |
|---|---|
| `--pos lat,lon,h` | approximate position for a warm start (km-level is enough) |
| `--cold` | cold start: no position needed |
| `--download` / `--nav FILE` | broadcast ephemeris from IGS (needs network) / a local RINEX 3 nav file |
| `--glonass` | also track GLONASS |
| `--coarse lat,lon` | cold start: rough position (hundreds of km) to search the likely visible satellites first; default: the system timezone's city |
| `--pty LINK` / `--nmea PATH` | NMEA (GGA, RMC, GSA) on a pty for gpsd / to a file or FIFO |
| `--secs N` | stop after N seconds (default: run until stopped) |
| `--log FILE` | log every chip message and solver event (large; useful for bug reports) |

## How it works

| file | what |
|---|---|
| `stpgps.py` | framing, escaping, checksum and payload scrambling of the `/dev/stpgps` protocol |
| `fwextract.py` | cuts the GPS-core firmware out of `libmnl.so` (SHA-256 checked) |
| `gpscfg.py` | the FE08 core configuration, built from named fields |
| `fwload.py` | bring-up: configuration, firmware download, core start; `--listen` decodes output |
| `ogps.py` | the engine: per-second aiding and tracking commands, clock and crystal model, cold start, NMEA |
| `pvt.py` | measurement decoding and solvers (least squares, EKF, cold single-epoch and multi-epoch) |
| `ephem.py` | GPS and GLONASS orbits, RINEX 3 navigation files, Doppler and az/el prediction |
| `navmsg.py` | GPS LNAV subframes and GLONASS strings to ephemerides; time from the sky |
| `gpsmon.py` | decoder and logger for a stock `mnld` session (from strace) |
| `logreplay.py` | re-runs the engine offline on a recorded `--log` session (solver testing, bug reports) |
| `correlate.py`, `replay.py` | reverse-engineering tools: message fields vs. predictions; replay of captured commands |
| `coldtest.sh` | cold-start test: wipes the saved state and logs the time to the first fix |
| `research/` | the experiments that identified the configuration fields and required messages |

### Transport (`/dev/stpgps`)

```
frame   = AA F0 | len:u16le | type:u16le | payload | cksum:u16le | AA 0F
len     : low 12 bits = bytes of type+payload+cksum; top 4 bits = rolling sequence (chip->host)
cksum   : 16-bit sum of the unescaped bytes from len to the end of the payload
escape  : 0xAA -> DE DF, 0xDE -> DE E0
payload : host->chip plain; chip->host plain = ror3(byte ^ 0xCA) (boot ROM replies FE31/FE32 plain)
```

### Bring-up

`FE08` configuration (sub-command 0x1a, 106 bytes), answered by `FE31`; 44 × `FE0C` (u16 index +
1400 bytes of firmware), each answered by `FE32 index`; `FE05 32 00 00 01 00`, answered by
`FE42` when the core runs (about 400 ms in total).

**FE08 0x1a configuration** (`gpscfg.py`): 52 × u16 after the id. Found by zeroing and varying
fields (`research/cfgbisect.py`, `research/cfgsem.py`):
- u16[8]=0x8000, [9]=0x8004, [10]=123, [11]=44 must have exactly these values (any other value
  tried: the core does not start). They look like a format or version signature.
- u16[30] is a flag (co-clock; non-zero needed).
- u16[29] is an enumerated clock mode (3 and 5 track; 0, 1, 2, 4 do not), next to
  u16[28]=25999, which looks like the 26 MHz reference in kHz.
- All other fields are tuning: the stock values are the defaults; `--minimal` leaves them zero
  and still gets fixes.

### Chip to host

**FE29 measurement** (57 bytes, about 1 Hz per channel):

| bytes | field |
|---|---|
| 0 / 1 / 2 | channel / sv (GLONASS: frequency channel + 8) / system (01 GPS, 10 GLONASS) |
| 3 | C/N0, dB-Hz |
| 5, 6 | flags; state: 0x40 locked (code + carrier, bit/frame sync, TOW), 0x20 pull-in, 0 lost |
| 23-24, 25 | **transmit time**: Z-count (TOW / 1.5 s, low 16 bits); high byte of the chip count since the start of that Z interval (0..23) |
| 27-28, 29-30 | chip count low 16 bits; fraction / 65536. tx = Z·1.5 s + (hi·65536 + chips + frac) / 1.023 MHz, in SV time. GLONASS: time of day + chips / 0.511 MHz |
| 31-32 | Doppler, Hz (int16) |
| 35-38 | accumulated carrier phase, 2^-16 cycles (u32, decreasing at positive Doppler) |

The receive time is the chip epoch (1 Hz, 32-bit millisecond counter in FE2C and FE28). All
channels of an epoch share it, so pseudorange = c·(T_epoch − tx) + c·dt_sv, and T_epoch's offset
is solved for.

**FE2B GPS subframe** (46 bytes): [2] sv, [6:46] 10 words × 4 bytes. The 24 data bits of a word
are in bytes 0, 3, 2 (byte 1 is a status byte). Word 1 starts with the 0x8B preamble; the data
are already corrected for the D30* inversion. Parity is not forwarded, and the chip forwards
subframes with bit errors, so an ephemeris is accepted only after two identical decodes.

**FE3B GLONASS string** (14 bytes): [0] frequency channel + 8, [2] = [3] string number, [4:14]
five little-endian u16 words = string bits 80..1, MSB first (data 80..9, Hamming code 8..1).
Strings that fail the Hamming check are dropped; an ephemeris needs two identical decodes. Times
are Moscow time (UTC+3); FE29 carries the GLONASS time of day as seconds mod 65536.

**FE2A GPS subframe sync** (10 bytes, every 6 s): [0] number of channels, [1:5] chip epoch
counter (hi16, lo16), [5:7] bitmask of channels that delivered the subframe, [7:9] TOW low 16
bits, [9] subframe id. **FE3A**: the GLONASS-string equivalent.

**FE42** (1 Hz): [12:14] crystal (DCXO) temperature reading, the input of the FE10 model.
**FE2C** epoch (32-bit ms counter at [0:4], hi16:lo16), **FE28** (same counter), **FE30**
clock/status (not decoded).

### Host to chip, every epoch

`FE00` next epoch counter, `FE10` crystal model, `FE01` time, one `FE02` per satellite, `FE03`
satellite lists, and housekeeping (`FE08` sub-commands, `FE1C`, `FE05 29 00`). Leaving out
messages one at a time (`research/ablate.sh`): without FE10, FE01 or `FE03 01 00` nothing
locks; the rest can be omitted.

**FE02** (20 bytes): `sv, sys, 01, 00, Z:u16, hi:u8, 00, chips:u16, doppler:i16, win1:u16,
win2:u16, flags[4]`: the predicted transmit time of the satellite at the next epoch (same
format as FE29) and its Doppler. Acquisition windows 402/512 (flags `00 00 01 00`, then
`02 00 01 00`); once the chip clock is known, 42/46 (`02 01 01 00`).

**FE03**: `00 n (sv, sys)…` is the aided list (satellites with an FE02); `01 2e (sv, sys)…` is
the autonomous search list (the core searches these by itself).

**FE01** (time): bias, drift, epoch counter, Z-count; bytes 4-5 are the frequency uncertainty
(warm about 0x0f8d, cold 0x1bb4; with 0 the core assumes its frequency is exact).

**FE10** (crystal model): a flag (`01` trusted, `00` not trusted), then 3 × (T:u16, 0:u32,
f:i16): the oscillator frequency offset f (about 1 Hz at L1 per unit) at the chip's temperature
reading T and at T ± 800. ogps sends f = c0 + 0.475·(T − 33438) and learns c0 in closed loop from
the residual Doppler, rate-limited to 5 Hz/s (larger steps make the carrier loops drop).

### Positioning (`pvt.py`)

- `EKF`: position, velocity, clock bias and drift, from pseudoranges and Doppler. Troposphere
  model; Klobuchar ionosphere with coefficients decoded from subframe 4 page 18; elevation and
  C/N0 weighting; innovation gating.
- `solve`: least-squares fix from code phase (mod 1 ms) near a known position.
- `solve_abs_robust`: cold single-epoch fix from the Earth's centre with 4+ satellites;
  20 ms navigation-bit ambiguities enumerated.
- `solve_abs_alt`: 3 satellites + Doppler + an altitude constraint.
- `solve_batch2`: several minutes of a static receiver with one continuous chip clock, for when
  fewer than 4 satellites are visible at once; needs 3+ satellites in total (2 can converge to
  a consistent but wrong place).

Transmit times are kept in extended precision (`numpy.longdouble`): absolute GPS seconds in
float64 resolve only 0.24 µs, i.e. 72 m of range. On platforms where `longdouble` is plain
double (32-bit ARM), fixes are correspondingly noisier.

### Almanac

The GPS almanac (subframe 4 pages 2-5/7-10, subframe 5 pages 1-24) is decoded from the sky and
cached in `almanac.json` next to the ephemeris cache (usable for weeks; pages are accepted when
the same content has been decoded twice; satellites can carry different uploads). It is only
used to aim the receiver, never for positioning:
- position known but no current orbits (offline, cache older than 4 h): FE02 aiding from the
  almanac (1-3 km orbit error, well inside the acquisition window). Outdoors, 5 satellites
  locked within 13 s and the first fix came after about 100 s (orbits decoded from the sky).
- cold start: with a coarse position (`--coarse`, or the timezone), the satellites likely above
  the horizon are searched first instead of in a fixed order.

### Cold start

Done as the stock engine does it, after capturing it with its aiding files deleted: no FE02 at
all; the autonomous search list in FE03, the larger frequency uncertainty in FE01, and an
untrusted FE10. Commanding blind searches with FE02 instead does not work: more than 8 of them
make the firmware assert (a whole-chip reset).

GPS time comes from the satellites (HOW TOW and the subframe-1 week) once two satellites agree.
Orbits come from subframes 1-3 after two identical decodes; single decodes are used only inside
the robust multi-epoch solver, which drops satellites that do not fit. A cold fix must be
confirmed by a second consistent solution before ogps switches to the aided loop.

## Known limits

- Only one device and firmware tested.
- FE30 and most FE08 sub-commands are not decoded; their stock values are sent.
- The firmware does not tolerate time-inconsistent aiding: replaying captured stock commands
  (`replay.py`) makes it assert.

## Tests

```sh
python3 -m pytest tests
```

The tests use a synthetic constellation and receiver, not recorded data: the transport codec,
GPS subframe and GLONASS string decoding, and the solvers (aided, cold single-epoch, bit
ambiguities, 3 satellites with altitude, multi-epoch batch).

## Contributing

Reports from other MT67xx devices are the most useful contribution. With the stock Android
stack running, capture its traffic and decode it:

```sh
strace -f -tt -y -xx -s 65535 -e trace=read,write -e signal=none -p $(pidof mnld) -o mnld.strace
python3 opengps/gpsmon.py -s mnld.strace
```

For problems with ogps itself, run it with `--log FILE`: the log can be re-run offline through
the engine with `opengps/logreplay.py FILE`, which reproduces the solver behaviour exactly.

Captures and logs include your position and time. Check what you share, or send only the summary
(`-s`) and the FE08/FE0C bring-up messages.

## Legal

The protocol was worked out from the device's own traffic, for interoperability. This
repository contains no MediaTek code and no firmware; the firmware has to be extracted from the
user's own device. MediaTek and other names are trademarks of their owners; this project is not
affiliated with them.

## Credits and references

- IS-GPS-200 (GPS LNAV message) and the GLONASS ICD (navigation strings).
- The GLONASS Hamming check masks C1..C7 are from RTKLIB (T. Takasu, BSD 2-clause); C_sigma follows
  the ICD (all 85 bits).
- Broadcast ephemerides: IGS, via the BKG GNSS data centre.
- [mt6572-connsys-modules](https://github.com/buddyjojo/mt6572-connsys-modules): the
  open connsys kernel driver.
