from datetime import datetime, timezone

from numpy import sqrt, pi, arccos, cos, sin
import numpy as np
import scipy.linalg as lin
from scipy.optimize import root_scalar
from matplotlib import pyplot as plt
from blocksim.satellite.Satellite import ASatellite, CircleSatellite
from blocksim.utils import itrf_to_teme, rotation_matrix, llavpa_to_itrf, rad, deg, itrf_to_azeld
from blocksim.constants import mu, Req


def mat_mvt_teme(t, tsync, pv0):
    t_epoch = (tsync - ASatellite.getInitialEpoch()).total_seconds()
    pv_teme = itrf_to_teme(t_epoch=t_epoch, pv_itrf=pv0)
    pos = pv_teme[:3]
    vel = pv_teme[3:]
    a = lin.norm(pos)
    sat_puls = sqrt(mu / a**3)

    n = np.cross(pos, vel)
    n /= lin.norm(n)

    R = rotation_matrix(axis=n, angle=sat_puls * t)

    return R


def mat_teme_itrf(t, tsync):
    """If R = mat_teme_itrf(t), and X is a position in TEME,
    then Y = R @ X  is the position in ITRF

    """
    T0 = 2451545.0
    DAY_S = 86400.0

    # Computing jd and fraction jd
    t_epoch = (tsync - ASatellite.getInitialEpoch()).total_seconds() + t
    epoch = t_epoch / 86400

    whole, fraction = divmod(epoch, 1.0)
    whole_jd = whole + 2433281.5

    jd_ut1 = whole_jd
    fraction_ut1 = fraction

    # Computing the rotation angle and axis
    t_teme = (jd_ut1 - T0 + fraction_ut1) / 36525.0
    g = 67310.54841 + (8640184.812866 + (0.093104 + (-6.2e-6) * t_teme) * t_teme) * t_teme
    theta = (jd_ut1 % 1.0 + fraction_ut1 + g / DAY_S % 1.0) % 1.0 * 2 * pi

    axe = np.array((0.0, 0.0, 1.0))

    # Computing the rotation matrix
    R = rotation_matrix(axis=axe, angle=-theta)

    return R


def main():
    h = 600e3
    r = Req + h

    tsync = datetime(2023, 1, 7, 13, 0, 0, tzinfo=timezone.utc)
    sat = CircleSatellite.fromOrbitalElements(
        name="sat", tsync=tsync, a=r, inc=45 * pi / 180, argp=0.5
    )
    Torb = sat.orbit_period.total_seconds()
    pv0 = sat.getGeocentricITRFPositionAt(0)
    M0 = pv0[:3] / r  # Normalize the satellite position

    lat = rad(43.60510103575826)
    lon = rad(1.4439216490854043)
    alt = 0
    elev = rad(20)
    pv_rx = llavpa_to_itrf(np.array([lon, lat, alt, 0, 0, 0]))
    pos_rx = pv_rx[:3] / lin.norm(pv_rx[:3])

    d = -Req * sin(elev) + sqrt(r**2 - Req**2 * cos(elev) ** 2)
    beta = arccos((d**2 + r**2 - Req**2) / (2 * r * d))
    s = sin(elev + beta)

    def fun(t):
        M1 = mat_teme_itrf(t, tsync) @ mat_mvt_teme(t, tsync, pv0) @ mat_teme_itrf(0, tsync).T @ M0

        # pv = sat.getGeocentricITRFPositionAt(t)
        # M1_ref = pv[:3]
        # err = M1 - M1_ref
        # assert lin.norm(err)<1e-8

        x = pos_rx.T @ M1
        return x - s

    res = root_scalar(f=fun, method="secant", x0=0, x1=Torb / 2)
    print(res)

    def check_criteria():
        tps = np.linspace(0, 5 * Torb, 501)
        elev = np.empty_like(tps)
        funval = np.empty_like(tps)
        for k in range(len(tps)):
            t = tps[k]
            pv = sat.getGeocentricITRFPositionAt(t)
            _, el, _, _, _, _ = itrf_to_azeld(pv_rx, pv)
            elev[k] = el
            funval[k] = fun(t)

        fig = plt.figure()
        axe = fig.add_subplot(211)
        axe.grid(True)
        axe.plot(tps, deg(elev))

        axe = fig.add_subplot(212, sharex=axe)
        axe.grid(True)
        axe.plot(tps, funval)

        plt.show()

    # check_criteria()


if __name__ == "__main__":
    main()
