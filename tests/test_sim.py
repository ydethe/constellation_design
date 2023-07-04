from numpy import pi

from constellation_design.simulate import simulate
from constellation_design.compute import nb_sat_min


def test_sgp4():
    inc = 37.5 * pi / 180
    nsat = 20
    npla = 20
    pha = 193
    alt = 607.4
    t_blind, nsat_max, tps_max, elev_mask = simulate(80 * pi / 180, inc, nsat, npla, pha, alt)
    print(f"t_blind: {t_blind:.1f} s")
    print(f"nsat_max: {nsat_max}")
    print(f"elev_mask: {elev_mask*180/pi:.2f} deg")


def test_optim():
    nb_sat_min()


if __name__ == "__main__":
    # test_sgp4()
    test_optim()
