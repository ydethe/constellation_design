from datetime import datetime, timezone

from numpy import sqrt, pi, arccos, cos, sin
import numpy as np
import scipy.linalg as lin
from matplotlib import pyplot as plt
from blocksim.satellite.Satellite import ASatellite, CircleSatellite
from blocksim.utils import itrf_to_teme, llavpa_to_itrf, rad, deg, itrf_to_azeld, rotation_matrix
from blocksim.constants import mu, Req


def get_axe_angle(R, prev_axe=None):
    # Angle determination
    s = np.trace(R)
    cos_angle = (s - 1) / 2
    if cos_angle > 1:
        cos_angle = 1
    if cos_angle < -1:
        cos_angle = -1
    angle = arccos(cos_angle)

    # Axe determination
    # https://math.stackexchange.com/questions/2074316/calculating-rotation-axis-from-rotation-matrix
    M = R + R.T + (1 - np.trace(R)) * np.eye(3)
    nc = lin.norm(M, axis=0)
    ir = np.argmax(nc)
    axe = M[ir, :]
    axe /= lin.norm(axe)

    # Matrix check
    P = rotation_matrix(angle, axe)
    if lin.norm(P - R) > 1e-3:
        angle = -angle
    P = rotation_matrix(angle, axe)

    # Continuity check
    if prev_axe is not None:
        if prev_axe @ axe < 0:
            axe = -axe
            angle = -angle

    return axe, angle


def mat_mvt_teme(t, tsync, pv0):
    t_epoch = (tsync - ASatellite.getInitialEpoch()).total_seconds()
    pv_teme = itrf_to_teme(t_epoch=t_epoch, pv_itrf=pv0)
    pos = pv_teme[:3]
    vel = pv_teme[3:]
    a = lin.norm(pos)
    sat_puls = sqrt(mu / a**3)
    angle = sat_puls * t

    n = np.cross(pos, vel)
    axe = n / lin.norm(n)
    kx, ky, kz = axe
    K = np.array([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]])
    K2 = np.einsum("i,j", axe, axe)

    if hasattr(t, "__iter__"):
        R = np.einsum("k,ij", cos(angle), np.eye(3))
        R += np.einsum("k,ij", sin(angle), K)
        R += np.einsum("k,ij", 1 - cos(angle), K2)
    else:
        R = cos(angle) * np.eye(3) + sin(angle) * K + (1 - cos(angle)) * K2

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
    angle = -(jd_ut1 % 1.0 + fraction_ut1 + g / DAY_S % 1.0) % 1.0 * 2 * pi

    axe = np.array([0.0, 0.0, 1.0])
    kx, ky, kz = axe
    K = np.array([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]])
    K2 = np.einsum("i,j", axe, axe)

    # Computing the rotation matrix
    if hasattr(t, "__iter__"):
        R = np.einsum("k,ij", cos(angle), np.eye(3))
        R += np.einsum("k,ij", sin(angle), K)
        R += np.einsum("k,ij", 1 - cos(angle), K2)
    else:
        R = cos(angle) * np.eye(3) + sin(angle) * K + (1 - cos(angle)) * K2

    return R


def main():
    h = 600e3
    r = Req + h

    tsync = datetime(2023, 1, 7, 13, 0, 0, tzinfo=timezone.utc)
    sat = CircleSatellite.fromOrbitalElements(
        name="sat", tsync=tsync, a=r, inc=45 * pi / 180, argp=0.5
    )
    sat.orbit_period.total_seconds()
    pv0 = sat.getGeocentricITRFPositionAt(0)
    M0 = mat_teme_itrf(0, tsync).T @ pv0[:3] / r  # Normalize the satellite position in TEME frame

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
        if hasattr(t, "__iter__"):
            R1 = mat_mvt_teme(t, tsync, pv0)
            R2 = mat_teme_itrf(t, tsync)
            Mt = np.einsum("ipj,p->ij", R1, M0)
            M1 = np.einsum("ip...,p...->i...", R2, Mt)
        else:
            M1 = mat_teme_itrf(t, tsync) @ mat_mvt_teme(t, tsync, pv0) @ M0

        x = pos_rx.T @ M1

        return x - s

    def check_criteria():
        tps = np.linspace(0, 86400, 2048)
        elev = np.empty_like(tps)
        funval2 = np.empty_like(tps)
        funval = fun(tps)
        for k in range(len(tps)):
            t = tps[k]
            pv = sat.getGeocentricITRFPositionAt(t)
            _, el, _, _, _, _ = itrf_to_azeld(pv_rx, pv)
            elev[k] = el
            funval2[k] = fun(t)

        fig = plt.figure()
        axe = fig.add_subplot(211)
        axe.grid(True)
        axe.plot(tps, deg(elev))

        axe = fig.add_subplot(212, sharex=axe)
        axe.grid(True)
        axe.plot(tps, funval, label="vector")
        axe.scatter(tps, funval2, label="point", marker=".")
        axe.legend()
        # axe.plot(tps, funval - funval2)

        plt.show()

    check_criteria()


if __name__ == "__main__":
    main()
