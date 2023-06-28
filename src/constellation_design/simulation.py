from datetime import datetime, timezone

import numpy as np
from blocksim.Simulation import Simulation
from blocksim.gnss.GNSSTracker import GNSSTracker
from blocksim.control.Route import Group
from blocksim.gnss.GNSSReceiver import GNSSReceiver
from blocksim.utils import rad
from blocksim.loggers.Logger import Logger
from blocksim.graphics.BFigure import FigureFactory
from blocksim.graphics import showFigures
from blocksim.graphics.GraphicSpec import AxeProjection
from blocksim.constants import Req


from .WalkerConstellation import WalkerConstellation


def build_sim(sma: float, inc: float, firstraan: float, t: int, p: int, f: int) -> Logger:
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

    tps = np.linspace(0, 3600, 30)
    sim.simulate(tps, progress_bar=True)

    log = sim.getLogger()
    print(log.getRawValue("tkr_vissat_n"))

    kmax = -1
    emax = -90
    for k in range(100):
        ec = np.max(log.getRawValue(f"tkr_obscoord_elev{k}"))
        if ec > emax:
            emax = ec
            kmax = k

    fig = FigureFactory.create()
    gs = fig.add_gridspec(1, 1)
    axe = fig.add_baxe(title="", spec=gs[0, 0], projection=AxeProjection.PLATECARREE)

    pt = (rec.lon, rec.lat)
    axe.plotDeviceReach(coord=pt, elev_min=tkr.elev_mask, sat_alt=sma - Req, color="blue")

    print(f"Highest sat: {satellites[kmax].getName() }")
    for s in satellites:
        traj = s.getTrajectoryFromLogger(log)
        axe.plot(traj, linewidth=4)

    showFigures()

    return log
