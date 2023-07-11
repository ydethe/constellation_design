from datetime import timedelta

import numpy as np
from numpy import pi
from pymoo.core.mixed import MixedVariableGA
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.variable import Real, Integer

from .simulate import simulate


# https://pymoo.org/customization/mixed.html


class NbSatMinimization(ElementwiseProblem):
    def __init__(self, **kwargs):
        self.alt_max_km = 1500
        self.max_lat_rad = 80 * pi / 180
        vars = {
            "nsatp": Integer(bounds=(1, 36)),
            "inc": Real(bounds=(60 * pi / 180, pi / 2)),
            "npla": Integer(bounds=(1, 36)),
            "pha": Integer(bounds=(0, 360)),
            "alt": Real(bounds=(500, self.alt_max_km)),
        }
        super().__init__(vars=vars, n_obj=1, n_ieq_constr=1, **kwargs)

    def _evaluate(self, X: dict, out: dict, *args, **kwargs):
        nsatp, inc, npla, pha, alt = X["nsatp"], X["inc"], X["npla"], X["pha"], X["alt"]
        nsat = nsatp * npla

        t_blind = np.empty(3)
        t_blind[0] = simulate(0, inc, nsat, npla, pha, alt)
        t_blind[1] = simulate(pi / 4, inc, nsat, npla, pha, alt)
        t_blind[2] = simulate(self.max_lat_rad, inc, nsat, npla, pha, alt)
        tb = timedelta(seconds=np.max(t_blind))
        print(
            f"[MIN]inc: {inc * 180 / pi:.1f}, nsat: {nsat}, npla: {npla},"
            f" pha: {pha}, alt: {alt :.1f},"
            f" t_blind: {tb}"
        )

        out["F"] = nsat
        out["G"] = t_blind @ t_blind - 3600**2


def nb_sat_min():
    problem = NbSatMinimization()

    algorithm = MixedVariableGA()

    res = minimize(
        problem,
        algorithm,
        termination=("n_evals", 10000),
        seed=1,
        verbose=True,
    )

    print("[MIN]Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
