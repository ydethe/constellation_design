from datetime import timedelta

import numpy as np
from numpy import pi

from pymoo.core.mixed import MixedVariableGA
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.variable import Real, Integer

from .simulate import simulate


# https://pymoo.org/customization/mixed.html


class MixedVariableProblem(ElementwiseProblem):
    def __init__(self, **kwargs):
        vars = {
            "inc": Real(bounds=(0, pi / 2)),
            # "nsat": Integer(bounds=(1, 1000)),
            "npla": Integer(bounds=(1, 36)),
            "pha": Integer(bounds=(0, 360)),
            "alt": Real(bounds=(500, 1500)),
        }
        self.nsat_per_plane = 1
        super().__init__(vars=vars, n_obj=1, n_ieq_constr=1, **kwargs)

    def _evaluate(self, X: dict, out: dict, *args, **kwargs):
        inc, npla, pha, alt = X["inc"], X["npla"], X["pha"], X["alt"]
        nsat = self.nsat_per_plane * npla

        t_blind = np.empty(3)
        t_blind[0], _, _, elev_mask = simulate(0, inc, nsat, npla, pha, alt)
        t_blind[1], _, _, _ = simulate(pi / 4, inc, nsat, npla, pha, alt)
        t_blind[2], _, _, _ = simulate(80 * pi / 180, inc, nsat, npla, pha, alt)
        tb = timedelta(seconds=np.max(t_blind))
        print(
            f"inc: {inc * 180 / pi:.1f}, npla: {npla},"
            f" pha: {pha}, alt: {alt :.1f}, elev: { elev_mask * 180 / pi:.1f},"
            f" t_blind: {tb}"
        )

        out["F"] = nsat
        out["G"] = t_blind @ t_blind - 3600**2


def optimize():
    problem = MixedVariableProblem()

    algorithm = MixedVariableGA(pop=10)

    res = minimize(
        problem,
        algorithm,
        termination=("n_evals", 10000),
        seed=1,
        verbose=True,
    )

    print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
