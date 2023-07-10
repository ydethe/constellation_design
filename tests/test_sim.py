from numpy import pi

from constellation_design.simulate import simulate
from constellation_design.compute import nb_sat_min


def test_sgp4():
    inc = 80 * pi / 180
    nsat = 1
    npla = 1
    pha = 193
    alt = 607.4
    t_blind = simulate(80 * pi / 180, inc, nsat, npla, pha, alt)
    print(f"t_blind: {t_blind:.1f} s")


def test_optim():
    nb_sat_min()


if __name__ == "__main__":
    # test_sgp4()
    test_optim()
