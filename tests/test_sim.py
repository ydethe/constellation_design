from datetime import timedelta

from numpy import pi, sqrt
import numpy as np
from scipy.optimize import minimize

from constellation_design.simulate import simulate

# from constellation_design.compute import nb_sat_min


def fun(X, nsat, npla, pha):
    inc, alt = X
    alt *= 1000
    t_blind = np.empty(3)
    t_blind[0] = simulate(0, inc, nsat, npla, pha, alt)
    t_blind[1] = simulate(pi / 4, inc, nsat, npla, pha, alt)
    t_blind[2] = simulate(80 * pi / 180, inc, nsat, npla, pha, alt)
    J = t_blind @ t_blind
    return J


def test_sgp4():
    inc = 80 * pi / 180
    nsat = 4
    npla = 4
    pha = 288
    alt = 894.8756129235549
    t_blind = np.empty(3)
    t_blind[0] = simulate(0, inc, nsat, npla, pha, alt)
    t_blind[1] = simulate(pi / 4, inc, nsat, npla, pha, alt)
    t_blind[2] = simulate(80 * pi / 180, inc, nsat, npla, pha, alt)
    print(f"t_blind: {np.max(t_blind):.1f} s")


def test_optim():
    global_min = None
    for nsatp in range(1, 36):
        for npla in range(1, 36):
            for pha in range(0, 360):
                nsat = nsatp * npla
                res = minimize(
                    fun=fun,
                    method="Powell",
                    x0=(75 * pi / 180, 1),
                    args=(nsat, npla, pha),
                    bounds=((60 * pi / 180, pi / 2), (0.60, 1.200)),
                )
                if global_min is None or res.fun < global_min:
                    global_min = res.fun
                    print(
                        nsat,
                        npla,
                        pha,
                        res.x[0] * 180 / pi,
                        res.x[1] * 1000,
                        timedelta(seconds=sqrt(res.fun)),
                    )

    # nb_sat_min()


if __name__ == "__main__":
    test_sgp4()
    # test_optim()
