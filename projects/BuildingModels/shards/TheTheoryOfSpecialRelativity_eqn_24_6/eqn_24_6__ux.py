from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_6__ux(self, c: float, u: float, v: float, **kwargs):
    # [.pyeqn] ux = u + v / (1 + (v * u) / c ** 2)
    result = []
    ux = (c**2 * u + c**2 * v + u**2 * v) / (c**2 + u * v)
    result.append(ux)
    return result
