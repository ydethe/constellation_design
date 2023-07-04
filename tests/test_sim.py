from constellation_design.simulate import simulate
from constellation_design.compute import optimize


def test_sgp4():
    t_blind, nsat_max, tps_max = simulate()
    print(t_blind, nsat_max, tps_max)


def test_optim():
    optimize()


if __name__ == "__main__":
    # test_sgp4()
    test_optim()
