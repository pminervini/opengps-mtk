# SPDX-License-Identifier: GPL-2.0-or-later
"""Position solvers on simulated FE29 measurements (synthetic constellation, Greenwich)."""
import os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim  # noqa: E402
import ephem, pvt  # noqa: E402

RX = np.array(sim.RX)


def test_fe29_transmit_time_roundtrip():
    nav = sim.constellation()
    ms, _ = sim.measurements(nav, sim.T0)
    assert len(ms) >= 6
    for m in ms:
        p = ephem.predict(nav, "G", m["sv"], sim.T0, sim.RX)
        assert abs(float(pvt.tx_from_meas(m, sim.T0) - sim.tx_time(p, sim.T0, tropo=True))) < 2e-9


def test_aided_solve():
    nav = sim.constellation()
    ms, t_rx = sim.measurements(nav, sim.T0)
    x0 = RX + np.array([20e3, -15e3, 5e3])                  # approximate position, 25 km off
    s = pvt.solve(ms, nav, t_rx, tuple(x0), min_el=0)
    assert np.linalg.norm(s["pos"] - RX) < 20                # troposphere model only on our side


def test_cold_single_epoch_from_earth_centre():
    nav = sim.constellation()
    ms, t_rx = sim.measurements(nav, sim.T0)
    s = pvt.solve_abs_robust(ms, nav, t_rx)
    assert np.linalg.norm(s["pos"] - RX) < 50


def test_cold_bit_ambiguity():
    nav = sim.constellation()
    ms, t_rx = sim.measurements(nav, sim.T0)
    for m in ms[:2]:                                        # two channels 9 navigation bits off
        m["tx_override"] = pvt.tx_from_meas(m, t_rx) - 0.180
    s = pvt.solve_abs_robust(ms, nav, t_rx)
    assert np.linalg.norm(s["pos"] - RX) < 50


def test_cold_three_satellites_altitude_doppler():
    nav = sim.constellation()
    ms, t_rx = sim.measurements(nav, sim.T0)
    ms = sorted(ms, key=lambda m: -ephem.predict(nav, "G", m["sv"], sim.T0, sim.RX)["el"])[:3]
    s = pvt.solve_abs_alt(ms, nav, t_rx, h0=70.0)
    assert np.linalg.norm(s["pos"] - RX) < 2000


def test_multi_epoch_batch_three_satellites():
    nav = sim.constellation()
    eps = []
    ms0, _ = sim.measurements(nav, sim.T0)
    sats = [m["sv"] for m in sorted(ms0, key=lambda m: -ephem.predict(nav, "G", m["sv"], sim.T0, sim.RX)["el"])[:3]]
    for k in range(180):                                    # 3 minutes, chip clock with a drift
        t = sim.T0 + k
        ms, _ = sim.measurements(nav, t, clock_bias=1e-4 + 2e-7 * k, sats=sats)
        eps.append((1000 * k, ms))
    s = pvt.solve_batch2(eps, nav, h0=70.0, step=6)
    assert s is not None and np.linalg.norm(s["pos"] - RX) < 3000
