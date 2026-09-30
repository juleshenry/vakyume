from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_12__v(
    self, L: float, d: float, delta_P: float, f: float, g: float, rho: float, **kwargs
):
    # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
    result = []
    v = -0.681202703290172 * sqrt(d * delta_P * g / (L * f * rho))
    result.append(v)
    v = 0.681202703290172 * sqrt(d * delta_P * g / (L * f * rho))
    result.append(v)
    return result
