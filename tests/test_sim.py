from constellation_design.simulation import build_sim


def test_sim():
    build_sim(
        sma=7000000.0,
        inc=85.0,
        firstraan=10.0,
        t=100,
        p=10,
        f=1,
    )

    # print(log.getRawValue("tkr_vissat_n"))


if __name__ == "__main__":
    test_sim()
