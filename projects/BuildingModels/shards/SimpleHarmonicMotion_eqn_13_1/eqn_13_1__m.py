from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_13_1__m(self, E: float, v: float, x: float, **kwargs):
    # [.pyeqn] E = 0.5 ** x ** 2 + 0.5 * m * v ** 2
    result = []
    m = 2.0 * E / v**2 - 2.0 / (2.0 ** (x**2) * v**2)
    result.append(m)
    return result
