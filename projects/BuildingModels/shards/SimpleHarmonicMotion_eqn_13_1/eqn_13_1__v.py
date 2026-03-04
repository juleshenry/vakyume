from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_13_1__v(self, E: float, m: float, x: float, **kwargs):
    # [.pyeqn] E = 0.5 ** x ** 2 + 0.5 * m * v ** 2
    result = []
    v = -1.4142135623731 * sqrt((E - 1 / 2.0 ** (x**2)) / m)
    result.append(v)
    v = 1.4142135623731 * sqrt((E - 1 / 2.0 ** (x**2)) / m)
    result.append(v)
    return result
