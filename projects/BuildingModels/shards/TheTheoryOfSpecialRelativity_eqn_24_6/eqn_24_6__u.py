from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_6__u(self, c: float, ux: float, v: float, **kwargs):
    # [.pyeqn] ux = u + v / (1 + (v * u) / c ** 2)
    result = []
    u = (
        -(c**2)
        + ux * v
        - sqrt((c**2 - 2 * c * v + ux * v) * (c**2 + 2 * c * v + ux * v))
    ) / (2 * v)
    result.append(u)
    u = (
        -(c**2)
        + ux * v
        + sqrt((c**2 - 2 * c * v + ux * v) * (c**2 + 2 * c * v + ux * v))
    ) / (2 * v)
    result.append(u)
    return result
