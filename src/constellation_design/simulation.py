from datetime import datetime, timezone

import more_itertools as mit
import numpy as np
from blocksim.Simulation import Simulation
from blocksim.gnss.GNSSTracker import GNSSTracker
from blocksim.control.Route import Group
from blocksim.gnss.GNSSReceiver import GNSSReceiver
from blocksim.utils import rad


from .WalkerConstellation import WalkerConstellation


def build_sim(sma: float, inc: float, firstraan: float, t: int, p: int, f: int) -> float:
    t0 = datetime(2023, 6, 27, 12, 0, 0, tzinfo=timezone.utc)
    sim = Simulation()

    rec = GNSSReceiver(
        name="rec",
        nsat=t,
        lat=rad(43.60510103575826),
        lon=rad(1.4439216490854043),
        alt=0,
        tsync=t0,
    )
    rec.algo = "no"

    const = WalkerConstellation("sim", sma, inc, firstraan, t, p, f)
    satellites = const.generate(tsync=t0)
    sim.addComputer(*satellites)

    tkr = GNSSTracker("tkr", t)
    tkr.elev_mask = rad(20)
    tkr.no_obs = True
    tkr.no_meas = True
    sim.addComputer(tkr)
    sim.addComputer(rec)

    # Then we do the connections
    nom_coord = ["px", "py", "pz", "vx", "vy", "vz"]

    grp_snames = []
    grp_inp = dict()
    for k, sat in enumerate(satellites):
        grp_inp["itrf%i" % k] = (6,)
        grp_snames.extend(["%s%i" % (n, k) for n in nom_coord])

    # The Group so defined let us gather all the outputs of the satellites
    # into one "wire" that feeds the tracker
    grp = Group(
        "grp",
        inputs=grp_inp,
        snames=grp_snames,
    )
    sim.addComputer(grp)

    for k, sat in enumerate(satellites):
        sim.connect("%s.itrf" % sat.getName(), "grp.itrf%i" % k)

    sim.connect("rec.realpos", "tkr.ueposition")
    sim.connect("grp.grouped", "tkr.state")
    sim.connect("tkr.measurement", "rec.measurements")
    sim.connect("tkr.ephemeris", "rec.ephemeris")

    tps = np.arange(0, 3 * sat.orbit_period.total_seconds(), 30)
    sim.simulate(tps, progress_bar=True)

    log = sim.getLogger()

    dt = tps[1] - tps[0]
    n = log.getRawValue("tkr_vissat_n")

    ind_nok = np.where(n == 0)[0]

    lg_max = -1
    g_max = [1]
    for group in mit.consecutive_groups(ind_nok):
        lgrp = list(group)
        lg = len(lgrp)
        if lg > lg_max:
            lg_max = lg
            g_max = lgrp.copy()

    worst_blind_time = (len(g_max) - 1) * dt

    return worst_blind_time
