from datetime import timedelta

import numpy as np
from numpy import pi, sin, cos, sqrt, arcsin
from blocksim.constants import Req
from pymoo.core.mixed import MixedVariableGA
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.variable import Real, Integer

from .simulate import simulate, compute_elevation_mask


# https://pymoo.org/customization/mixed.html


class NbSatMinimization(ElementwiseProblem):
    def __init__(self, **kwargs):
        self.alt_max_km = 1500
        self.max_lat_rad = 80 * pi / 180
        elev = compute_elevation_mask(self.alt_max_km)
        rad = Req + self.alt_max_km * 1e3
        dist = -Req * sin(elev) + sqrt(rad**2 - Req**2 * cos(elev) ** 2)

        a = arcsin(dist * cos(elev) / rad)
        inc_min_rad = self.max_lat_rad - a
        print(
            f"Min. inclination for {self.alt_max_km:.1f} km max altitude:"
            f" {inc_min_rad*180/pi:.2f} deg"
        )

        vars = {
            "inc": Real(bounds=(inc_min_rad, pi / 2)),
            "npla": Integer(bounds=(1, 36)),
            "pha": Integer(bounds=(0, 360)),
            "alt": Real(bounds=(500, self.alt_max_km)),
        }
        self.nsat_per_plane = 1
        super().__init__(vars=vars, n_obj=1, n_ieq_constr=1, **kwargs)

    def _evaluate(self, X: dict, out: dict, *args, **kwargs):
        inc, npla, pha, alt = X["inc"], X["npla"], X["pha"], X["alt"]
        nsat = self.nsat_per_plane * npla

        t_blind = np.empty(3)
        t_blind[0], _, _, elev_mask = simulate(0, inc, nsat, npla, pha, alt)
        t_blind[1], _, _, _ = simulate(pi / 4, inc, nsat, npla, pha, alt)
        t_blind[2], _, _, _ = simulate(self.max_lat_rad, inc, nsat, npla, pha, alt)
        tb = timedelta(seconds=np.max(t_blind))
        print(
            f"[MIN]inc: {inc * 180 / pi:.1f}, npla: {npla},"
            f" pha: {pha}, alt: {alt :.1f}, elev: { elev_mask * 180 / pi:.1f},"
            f" t_blind: {tb}"
        )

        out["F"] = nsat
        out["G"] = t_blind @ t_blind - 3600**2


def nb_sat_min():
    problem = NbSatMinimization()

    algorithm = MixedVariableGA(pop=10)

    res = minimize(
        problem,
        algorithm,
        termination=("n_evals", 10000),
        seed=1,
        verbose=True,
    )

    print("[MIN]Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
