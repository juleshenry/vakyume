from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_3__installed_costs(self, hp: float, **kwargs):
    # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
    result = []
    installed_costs = 13482.9087908759 * hp ** (9 / 20)
    result.append(installed_costs)
    return result
