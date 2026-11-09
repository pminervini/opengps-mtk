# SPDX-License-Identifier: GPL-2.0-or-later
"""Broadcast ephemeris (RINEX 3 navigation files) -> satellite position, velocity, clock.

GPS: IS-GPS-200 Kepler model. GLONASS: PZ-90 state vector integrated with RK4 (GLONASS ICD).
Times are GPS seconds since 1980-01-06 ("gpst") unless noted.
"""
import gzip, math, calendar

GM = 3.986005e14
OMEGA_E = 7.2921151467e-5
C = 299792458.0
F_REL = -4.442807633e-10
L1 = 1575.42e6
LEAP = 18                                   # GPS - UTC, seconds (valid since 2017)
GPS_EPOCH = calendar.timegm((1980, 1, 6, 0, 0, 0))
WGS_A, WGS_F = 6378137.0, 1 / 298.257223563
# GLONASS (PZ-90.11)
GL_GM, GL_AE, GL_J2, GL_OMEGA = 3.986004418e14, 6378136.0, 1.08262575e-3, 7.292115e-5


def unix_to_gpst(u):
    return u - GPS_EPOCH + LEAP


def gpst_to_unix(g):
    return g + GPS_EPOCH - LEAP


def _f(s):
    s = s.strip().replace("D", "E").replace("d", "e")
    return float(s) if s else 0.0


def _vals(line, start=4, n=4):
    return [_f(line[start + 19 * i:start + 19 * (i + 1)]) for i in range(n)]


def load_rinex(path):
    """-> {'G': {prn: [eph...]}, 'R': {slot: [eph...]}} (dicts of floats)."""
    op = gzip.open if path.endswith(".gz") else open
    lines = op(path, "rt").read().splitlines()
    i = next(k for k, l in enumerate(lines) if "END OF HEADER" in l) + 1
    out = {"G": {}, "R": {}}
    while i < len(lines):
        l = lines[i]
        sysc = l[:1]
        nlines = {"G": 7, "E": 7, "C": 7, "J": 7, "I": 7, "R": 3, "S": 3}.get(sysc, 0)
        if not nlines:
            i += 1
            continue
        if sysc == "R" and i + 4 < len(lines) and lines[i + 4][:1] == " ":
            nlines = 4                       # RINEX 3.05 GLONASS has a 4th line
        prn = int(l[1:3])
        y, mo, d, h, mi, s = int(l[4:8]), int(l[9:11]), int(l[12:14]), int(l[15:17]), int(l[18:20]), int(l[21:23])
        head = _vals(l, 23, 3)
        body = [v for k in range(1, nlines + 1) for v in _vals(lines[i + k])]
        i += nlines + 1
        t_unix = calendar.timegm((y, mo, d, h, mi, s, 0, 0, 0))
        if sysc == "G":
            e = dict(prn=prn, toc=t_unix - GPS_EPOCH, af0=head[0], af1=head[1], af2=head[2])
            keys = ["iode", "crs", "dn", "m0", "cuc", "e", "cus", "sqrta", "toe", "cic", "omega0", "cis",
                    "i0", "crc", "w", "omegadot", "idot", "l2c", "week", "l2p", "acc", "health", "tgd", "iodc",
                    "ttm", "fit"]
            e.update(zip(keys, body))
            e["toe_g"] = e["week"] * 604800 + e["toe"]
            out["G"].setdefault(prn, []).append(e)
        elif sysc == "R":
            # epoch in UTC(SU) ~ UTC; store as gpst
            e = dict(slot=prn, tb=unix_to_gpst(t_unix), taun=-head[0], gamman=head[1],
                     x=body[0] * 1e3, vx=body[1] * 1e3, ax=body[2] * 1e3, health=body[3],
                     y=body[4] * 1e3, vy=body[5] * 1e3, ay=body[6] * 1e3, k=int(body[7]),
                     z=body[8] * 1e3, vz=body[9] * 1e3, az=body[10] * 1e3)
            out["R"].setdefault(prn, []).append(e)
    return out


def pick(ephs, t, key):
    return min(ephs, key=lambda e: abs(t - e[key])) if ephs else None


def gps_sat(e, t):
    """GPS sat at transmit time t (gpst) -> (pos ECEF m, vel m/s, clock bias s, clock drift s/s)."""
    t = float(t)                        # may arrive in extended precision; mm-level here
    a = e["sqrta"] ** 2
    tk = t - e["toe_g"]
    n = math.sqrt(GM / a ** 3) + e["dn"]
    M = e["m0"] + n * tk
    E = M
    for _ in range(10):
        E = M + e["e"] * math.sin(E)
    sE, cE = math.sin(E), math.cos(E)
    Edot = n / (1 - e["e"] * cE)
    v = math.atan2(math.sqrt(1 - e["e"] ** 2) * sE, cE - e["e"])
    phi = v + e["w"]
    s2, c2 = math.sin(2 * phi), math.cos(2 * phi)
    du = e["cus"] * s2 + e["cuc"] * c2
    dr = e["crs"] * s2 + e["crc"] * c2
    di = e["cis"] * s2 + e["cic"] * c2
    u = phi + du
    r = a * (1 - e["e"] * cE) + dr
    i = e["i0"] + e["idot"] * tk + di
    xo, yo = r * math.cos(u), r * math.sin(u)
    Om = e["omega0"] + (e["omegadot"] - OMEGA_E) * tk - OMEGA_E * e["toe"]
    cO, sO, ci, si = math.cos(Om), math.sin(Om), math.cos(i), math.sin(i)
    pos = (xo * cO - yo * ci * sO, xo * sO + yo * ci * cO, yo * si)
    # velocity (analytic)
    vdot = Edot * math.sqrt(1 - e["e"] ** 2) / (1 - e["e"] * cE)
    udot = vdot * (1 + 2 * (e["cus"] * c2 - e["cuc"] * s2))
    rdot = a * e["e"] * sE * Edot + 2 * vdot * (e["crs"] * c2 - e["crc"] * s2)
    idot = e["idot"] + 2 * vdot * (e["cis"] * c2 - e["cic"] * s2)
    Omdot = e["omegadot"] - OMEGA_E
    xod = rdot * math.cos(u) - r * udot * math.sin(u)
    yod = rdot * math.sin(u) + r * udot * math.cos(u)
    vel = (xod * cO - yod * ci * sO + yo * si * sO * idot - pos[1] * Omdot,
           xod * sO + yod * ci * cO - yo * si * cO * idot + pos[0] * Omdot,
           yod * si + yo * ci * idot)
    dt = t - e["toc"]
    clk = e["af0"] + e["af1"] * dt + e["af2"] * dt * dt + F_REL * e["e"] * e["sqrta"] * sE - e["tgd"]
    drift = e["af1"] + 2 * e["af2"] * dt
    return pos, vel, clk, drift


def _gl_deriv(s, acc):
    x, y, z, vx, vy, vz = s
    r2 = x * x + y * y + z * z
    r = math.sqrt(r2)
    g = -GL_GM / r ** 3
    j = 1.5 * GL_J2 * GL_GM * GL_AE ** 2 / r ** 5
    zz = 5 * z * z / r2
    w2 = GL_OMEGA ** 2
    ax = g * x - j * x * (1 - zz) + w2 * x + 2 * GL_OMEGA * vy + acc[0]
    ay = g * y - j * y * (1 - zz) + w2 * y - 2 * GL_OMEGA * vx + acc[1]
    az = g * z - j * z * (3 - zz) + acc[2]
    return (vx, vy, vz, ax, ay, az)


def glo_sat(e, t):
    """GLONASS sat at transmit time t (gpst) -> (pos, vel, clock bias s, drift)."""
    t = float(t)
    s = (e["x"], e["y"], e["z"], e["vx"], e["vy"], e["vz"])
    acc = (e["ax"], e["ay"], e["az"])
    dt = t - e["tb"]
    h = 60.0 if dt >= 0 else -60.0
    while abs(dt) > 1e-9:
        st = h if abs(dt) > abs(h) else dt
        k1 = _gl_deriv(s, acc)
        k2 = _gl_deriv([s[i] + st / 2 * k1[i] for i in range(6)], acc)
        k3 = _gl_deriv([s[i] + st / 2 * k2[i] for i in range(6)], acc)
        k4 = _gl_deriv([s[i] + st * k3[i] for i in range(6)], acc)
        s = [s[i] + st / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(6)]
        dt -= st
    clk = -e["taun"] + e["gamman"] * (t - e["tb"])
    return tuple(s[:3]), tuple(s[3:]), clk, e["gamman"]


def geo2ecef(lat, lon, h):
    lat, lon = math.radians(lat), math.radians(lon)
    e2 = WGS_F * (2 - WGS_F)
    N = WGS_A / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    return ((N + h) * math.cos(lat) * math.cos(lon), (N + h) * math.cos(lat) * math.sin(lon),
            (N * (1 - e2) + h) * math.sin(lat))


def ecef2geo(x, y, z):
    e2 = WGS_F * (2 - WGS_F)
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1 - e2))
    for _ in range(6):
        N = WGS_A / math.sqrt(1 - e2 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - N
        lat = math.atan2(z, p * (1 - e2 * N / (N + h)))
    return math.degrees(lat), math.degrees(lon), h


def azel(rx, sat):
    lat, lon, _ = ecef2geo(*rx)
    la, lo = math.radians(lat), math.radians(lon)
    d = [sat[i] - rx[i] for i in range(3)]
    e = -math.sin(lo) * d[0] + math.cos(lo) * d[1]
    n = -math.sin(la) * math.cos(lo) * d[0] - math.sin(la) * math.sin(lo) * d[1] + math.cos(la) * d[2]
    u = math.cos(la) * math.cos(lo) * d[0] + math.cos(la) * math.sin(lo) * d[1] + math.sin(la) * d[2]
    return math.degrees(math.atan2(e, n)) % 360, math.degrees(math.atan2(u, math.hypot(e, n)))


def predict(nav, sysc, prn, t_rx, rx):
    """Receiver at rx (ECEF) receiving at t_rx (gpst): geometric range, range rate, sat clock, az/el.
    Handles light time and Earth rotation (Sagnac)."""
    t_rx = float(t_rx)
    if sysc == "G":
        e = pick(nav["G"].get(prn, []), t_rx, "toe_g")
        f = gps_sat
    else:
        e = pick(nav["R"].get(prn, []), t_rx, "tb")
        f = glo_sat
    if not e:
        return None
    tau = 0.075
    for _ in range(3):
        pos, vel, clk, drift = f(e, t_rx - tau)
        th = OMEGA_E * tau
        p = (pos[0] * math.cos(th) + pos[1] * math.sin(th), -pos[0] * math.sin(th) + pos[1] * math.cos(th), pos[2])
        rng = math.dist(p, rx)
        tau = rng / C
    los = [(p[i] - rx[i]) / rng for i in range(3)]
    rr = sum(los[i] * vel[i] for i in range(3))       # range rate, m/s (receiver static)
    az, el = azel(rx, p)
    return dict(range=rng, rate=rr, clk=clk, drift=drift, az=az, el=el, pos=p, vel=vel, eph=e)
