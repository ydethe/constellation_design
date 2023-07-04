from numpy import pi

from constellation_design.simulate import simulate
from constellation_design.compute import optimize


def test_sgp4():
    inc = 1.6 * pi / 180
    nsat = 33
    npla = 33
    pha = 212
    alt = 1241.0
    t_blind, nsat_max, tps_max, elev_mask = simulate(80 * pi / 180, inc, nsat, npla, pha, alt)
    print(f"t_blind: {t_blind:.1f} s")
    print(f"nsat_max: {nsat_max}")
    print(f"elev_mask: {elev_mask*180/pi:.2f} deg")


def test_optim():
    optimize()


if __name__ == "__main__":
    # test_sgp4()
    test_optim()
