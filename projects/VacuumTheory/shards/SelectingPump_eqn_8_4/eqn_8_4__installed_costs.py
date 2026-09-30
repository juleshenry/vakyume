from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_4__installed_costs(self, hp: float, **kwargs):
    # [.pyeqn] installed_costs = 26000 * (hp / 10) ** 0.4
    result = []
    installed_costs = 10350.7864343909 * hp ** (2 / 5)
    result.append(installed_costs)
    return result
