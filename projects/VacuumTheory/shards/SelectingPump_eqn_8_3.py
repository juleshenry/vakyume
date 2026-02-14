from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SelectingPump:
    @kwasak_static
    def eqn_8_3(hp=None, installed_costs=None, **kwargs):
        return

    @staticmethod
    def eqn_8_3__hp(installed_costs: float, **kwargs):
        # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
        # Error during Sympy solve: Sympy solve failed
        def func(hp):
            # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
            return eval("(38000 * (x / 10) ** 0.45) - (installed_costs)".replace('x', str(hp)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_8_3__installed_costs(hp: float, **kwargs):
        # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
        result = []
        installed_costs = 13482.9087908759*hp**(9/20)
        result.append(installed_costs)
        return result

