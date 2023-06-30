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
    assert lin.norm(P - R) < 1e-3

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
    kx, ky, kz = n / lin.norm(n)
    K = np.array([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]])
    K2 = K @ K

    if hasattr(t, "__iter__"):
        R = np.einsum("k,ij", np.ones_like(angle), np.eye(3))
        R += np.einsum("k,ij", sin(angle), K)
        R += np.einsum("k,ij", 1 - cos(angle), K2)

        for k in range(len(t)):
            err = lin.norm(R[:, :, k] @ R[:, :, k].T - np.eye(3))
            assert err < 1e-6
    else:
        R = np.eye(3) + sin(angle) * K + (1 - cos(angle)) * K2

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
    angle = (jd_ut1 % 1.0 + fraction_ut1 + g / DAY_S % 1.0) % 1.0 * 2 * pi

    kx, ky, kz = 0.0, 0.0, 1.0
    K = np.array([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]])
    K2 = K @ K

    # Computing the rotation matrix
    if hasattr(t, "__iter__"):
        R = np.einsum("k,ij", np.ones_like(angle), np.eye(3))
        R += np.einsum("k,ij", -sin(angle), K)
        R += np.einsum("k,ij", 1 - cos(angle), K2)

        for k in range(len(t)):
            err = lin.norm(R[:, :, k] @ R[:, :, k].T - np.eye(3))
            assert err < 1e-6
    else:
        R = np.eye(3) - sin(angle) * K + (1 - cos(angle)) * K2

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
            for k in range(len(t)):
                err = lin.norm(R1[:, :, k] @ M0 - Mt[:, k])
                assert err < 1e-6
            M1 = np.einsum("ip...,p...->i...", R2, Mt)
            for k in range(len(t)):
                err = lin.norm(R2[:, :, k] @ R1[:, :, k] @ M0 - M1[:, k])
                assert err < 1e-6
        else:
            M1 = mat_teme_itrf(t, tsync) @ mat_mvt_teme(t, tsync, pv0) @ M0

        # pv = sat.getGeocentricITRFPositionAt(t)
        # M1_ref = pv[:3]
        # err = M1 - M1_ref
        # assert lin.norm(err)<1e-8

        x = pos_rx.T @ M1
        if hasattr(t, "__iter__"):
            for k in range(len(t)):
                err = np.abs(x[k] - pos_rx.T @ R2[:, :, k] @ R1[:, :, k] @ M0)
                assert err < 1e-6

        return x - s

    def check_criteria():
        tps = np.linspace(0, 86400 * 5, 128 * 5 + 1)
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
        axe.plot(tps, funval2, label="point")
        axe.legend()
        # axe.plot(tps, funval - funval2)

        plt.show()

    check_criteria()


if __name__ == "__main__":
    main()
