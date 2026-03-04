from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_13_1__x(self, E: float, m: float, v: float, **kwargs):
    # [.pyeqn] E = 0.5 ** x ** 2 + 0.5 * m * v ** 2
    result = []
    x = -1.20112240878645 * sqrt(log(2.0 / (2.0 * E - m * v**2)))
    result.append(x)
    x = 1.20112240878645 * sqrt(log(2.0 / (2.0 * E - m * v**2)))
    result.append(x)
    return result
