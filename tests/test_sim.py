from constellation_design.simulation import build_sim


def test_sim():
    worst_blind_time = build_sim(
        sma=7000000.0,
        inc=85.0,
        firstraan=10.0,
        t=100,
        p=10,
        f=1,
    )

    print(worst_blind_time)


if __name__ == "__main__":
    test_sim()
