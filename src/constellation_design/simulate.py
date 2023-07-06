from datetime import datetime, timezone
from dataclasses import dataclass
from enum import Enum

import numpy as np
from numpy import arcsin, pi, log10
from numpy.polynomial import Polynomial
from blocksim.satellite.Satellite import (
    CircleSatellite,
    generateWalkerDeltaConstellation,
)
from blocksim.utils import llavpa_to_itrf, itrf_to_azeld
from blocksim.constants import Req, kb, c as clum


class EventType(Enum):
    RISE = "rise"
    SET = "set"


@dataclass
class Event:
    satellite: int
    type: EventType
    date: float


def analyse_timeline(init, events, total_sim_time):
    tl = []
    for ksat, sat in enumerate(events):
        for e in sat:
            tl.append(Event(satellite=ksat, type=EventType.RISE, date=e["rise"]))
            tl.append(Event(satellite=ksat, type=EventType.SET, date=e["set"]))
    tl.sort(key=lambda x: x.date)

    nsat = len(np.where(init > 0)[0])
    blind_e = Event(satellite=-1, type=EventType.RISE, date=0)
    t_blind = 0
    e: Event
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
        print(e,nsat,t_blind)

    if nsat_max == 0:
        t_blind = total_sim_time

    return t_blind, nsat_max


def compute_elevation_mask(alt_km):
    # Parametres du probleme
    # ======================
    px = 36  # dBm
    NF = 3  # dB
    z1_J1 = 1.616339347
    eta = 0.7
    wl = clum / (2e9)
    cn0_lim = 46.46595211
    alpha = 0.0305

    # Calculs des constantes K et Q
    # =============================
    sma = Req + alt_km * 1e3
    K = (
        -36
        + px
        - NF
        + 10
        * log10((eta**2 * wl**2 * z1_J1**4 * sma**4) / (16 * pi**2 * kb * 290 * Req**4))
    )
    Q = 10 ** ((K - cn0_lim) / 20)

    # Calcul de l'élévation
    # =====================
    p = Polynomial(
        [
            alpha,
            2,
            -(
                (2 * alpha**2 - 1) * sma**2
                + (1 - 2 * alpha**2) * Req**2
                + 2 * alpha * Q * Req
                + Q**2
            )
            / (alpha * sma**2 - alpha * Req**2),
            -(4 * alpha * sma**2 - 4 * alpha * Req**2 + 2 * Q * Req)
            / (alpha * sma**2 - alpha * Req**2),
            ((alpha**2 - 2) * sma**2 + (2 - alpha**2) * Req**2 + 2 * alpha * Q * Req)
            / (alpha * sma**2 - alpha * Req**2),
            (2 * alpha * sma**2 - 2 * alpha * Req**2 + 2 * Q * Req)
            / (alpha * sma**2 - alpha * Req**2),
            1 / alpha,
        ]
    )

    rts = p.roots()
    r0 = rts[np.where((rts > arcsin(alpha)) & (rts < 1) & (np.abs(np.imag(rts)) < 1e-9))[0]]
    assert len(r0) == 1
    elev_mask = arcsin(r0[0])

    return elev_mask


def simulate(lat, inc, nsat, npla, pha, alt_km):
    elev_mask = compute_elevation_mask(alt_km)

    t0 = datetime(2023, 6, 27, 12, 0, 0, tzinfo=timezone.utc)
    firstraan = 0.0
    lon = 0.0
    # tps_max = 5 * 86400
    tps_max = 20000
    sma = Req + alt_km * 1e3

    satellites = generateWalkerDeltaConstellation(
        "sim", sma, inc, firstraan, nsat, npla, pha, t0, prop=CircleSatellite
    )
    obs = llavpa_to_itrf((lon, lat, 0, 0, 0, 0))
    events = list()
    init = list()
    sat: CircleSatellite
    for sat in satellites:
        pv_sat = sat.getGeocentricITRFPositionAt(0)
        _, el0, _, _, _, _ = itrf_to_azeld(obs, pv_sat)
        init.append(el0 - elev_mask)
        events.append(
            sat.find_events(obs, t0=0, t1=tps_max, elevation=elev_mask),
        )

    t_blind, nsat_max = analyse_timeline(np.array(init), events, tps_max)

    return t_blind, nsat_max, tps_max, elev_mask
