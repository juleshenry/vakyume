from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_12__delta_P(
    self, L: float, d: float, f: float, g: float, rho: float, v: float, **kwargs
):
    # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
    result = []
    delta_P = 2.155 * L * f * rho * v**2 / (d * g)
    result.append(delta_P)
    return result
