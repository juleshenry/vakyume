from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SelectingPump:
    @kwasak_static
    def eqn_8_02(hp: float = None, installed_costs: float = None, **kwargs):
        return

    @staticmethod
    def eqn_8_02__hp(installed_costs: float):
        # [.pyeqn] installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        hp = 9.18273645546364e-9 * installed_costs**2
        result.append(hp)
        return result

    @staticmethod
    def eqn_8_02__installed_costs(hp: float):
        # [.pyeqn] installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        installed_costs = 10435.5162785557 * sqrt(hp)
        result.append(installed_costs)
        return result


