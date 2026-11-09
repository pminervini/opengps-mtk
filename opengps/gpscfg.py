# SPDX-License-Identifier: GPL-2.0-or-later
"""FE08 sub-command 0x1a: GPS-core configuration (106 bytes = id u16 + 52 x u16).

Built from named fields instead of a copied blob. Values are what stock mnld sends on the
MT6750 tablet (NVRAM: co-clock, 26 MHz). Measured on the device (tools/gps/opengps/README.md):
  * the core does not start (no FE42) if any of u16[8..11] or u16[30] is zero
  * with only those five, it starts but tracks nothing; u16[29] is also needed for tracking
  * everything else is tuning (zeroed: still starts, tracks and gets fixes)
"""
import struct

REQUIRED = {
    8: 0x8000,     # mode flags (needed to start)
    9: 0x8004,     # mode flags (needed to start)
    10: 123,       # needed to start; meaning unknown
    11: 44,        # needed to start; meaning unknown
    29: 3,         # needed for tracking; next to the clock frequency -> clock type (3 = co-clock DCXO?)
    30: 1,         # needed to start; co-clock flag (COMBO_IOC_CO_CLOCK_FLAG returns 1 on this board)
}

TUNING = {          # stock mnld defaults; optional
    3: 451, 4: 677, 7: 0x7800, 14: 5, 15: 200, 16: 0x3C6F, 17: 0x8349, 18: 1511,
    20: 210, 21: 183, 22: 40, 23: 19032, 24: 30, 25: 14274, 26: 40, 27: 19032,
    28: 25999,     # 26 MHz reference clock in kHz (-1)? (NVRAM: 26000000 Hz)
    31: 40, 32: 255, 33: 128, 34: 71, 35: 100, 36: 80, 37: 216, 38: 80, 39: 187, 40: 3,
    41: 60, 42: 111, 43: 137, 44: 136, 45: 2, 46: 251, 47: 1,
    49: 0x4948, 50: 0x4B4A, 51: 0x4D4C, 52: 0x4F4E,   # "HIJKLMNO" (test pattern / padding)
}


def build(minimal=False):
    f = [0] * 53
    f[0] = 0x001A                   # sub-command id
    f_upd = dict(REQUIRED)
    if not minimal:
        f_upd.update(TUNING)
    for k, v in f_upd.items():
        f[k] = v
    return struct.pack("<53H", *f)


if __name__ == "__main__":
    print(build().hex(" "))
