from datetime import datetime, timezone
from dataclasses import dataclass
from enum import Enum

import numpy as np
from numpy import pi, sqrt
import rich.progress as rp
from blocksim.satellite import (
    CircleSatellite,
    generateWalkerDeltaConstellation,
)
from blocksim.utils import rad, llavpa_to_itrf, itrf_to_azeld
from blocksim.constants import mu


class EventType(Enum):
    RISE = "rise"
    SET = "set"


@dataclass
class Event:
    satellite: int
    type: EventType
    date: float


def analyse_timeline(init, events):
    tl = []
    for ksat, sat in enumerate(events):
        for e in sat:
            tl.append(Event(satellite=ksat, type=EventType.RISE, date=e["rise"]))
            tl.append(Event(satellite=ksat, type=EventType.SET, date=e["set"]))
    tl.sort(key=lambda x: x.date)

    nsat = len(np.where(init > 0)[0])
    blind_e: EventType = None
    t_blind = 0
    e: EventType
    nsat_max = nsat
    for e in tl:
        if e.type == EventType.RISE:
            if nsat == 0:
                t_blind += e.date - blind_e.date
            nsat += 1
        else:
            nsat -= 1
            if nsat == 0:
                blind_e = e
            elif nsat < 0:
                raise AssertionError
        if nsat_max < nsat:
            nsat_max = nsat

    return t_blind, nsat_max


def test_sgp4():
    t0 = datetime(2023, 6, 27, 12, 0, 0, tzinfo=timezone.utc)
    sma = 7000000.0
    inc = 85.0
    firstraan = 10.0
    t = 100
    p = 10
    f = 1
    lat = rad(43.60510103575826)
    lon = rad(1.4439216490854043)
    elev_mask = rad(20)
    tps_max = 5 * 86400

    satellites = generateWalkerDeltaConstellation(
        "sim", sma, inc, firstraan, t, p, f, t0, prop=CircleSatellite
    )
    2 * pi * sqrt(sma**3 / mu)
    obs = llavpa_to_itrf((lon, lat, 0, 0, 0, 0))
    events = list()
    init = list()
    for sat in rp.track(satellites):
        pv_sat = sat.getGeocentricITRFPositionAt(0)
        _, el0, _, _, _, _ = itrf_to_azeld(obs, pv_sat)
        init.append(el0 - elev_mask)
        events.append(
            sat.find_events(obs, t0=0, t1=tps_max, elevation=elev_mask),
        )

    t_blind, nsat_max = analyse_timeline(np.array(init), events)
    print(t_blind, nsat_max, tps_max)


if __name__ == "__main__":
    test_sgp4()
