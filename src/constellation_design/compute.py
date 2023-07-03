from pymoo.core.mixed import MixedVariableGA
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.variable import Real, Integer


# https://pymoo.org/customization/mixed.html


class MixedVariableProblem(ElementwiseProblem):
    def __init__(self, **kwargs):
        vars = {
            "inc": Real(bounds=(0, 180)),
            "nsat": Integer(bounds=(1, 10000)),
            "npla": Integer(bounds=(1, 36)),
            "pha": Real(bounds=(0, "npla")),
            "alt": Real(bounds=(500, 1500)),
        }
        super().__init__(vars=vars, n_obj=1, **kwargs)

    def _evaluate(self, X: dict, out: dict, *args, **kwargs):
        inc, nsat, npla, pha, alt = X["inc"], X["nsat"], X["npla"], X["pha"], X["alt"]

        inc + nsat + npla + pha + alt

        out["F"] = nsat


problem = MixedVariableProblem()

algorithm = MixedVariableGA(pop=10)

res = minimize(problem, algorithm, termination=("n_evals", 10000), seed=1, verbose=False)

print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
