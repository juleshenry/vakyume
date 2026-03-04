from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_6__v(self, c: float, u: float, ux: float, **kwargs):
    # [.pyeqn] ux = u + v / (1 + (v * u) / c ** 2)
    result = []
    v = c**2 * (-u + ux) / (c**2 + u**2 - u * ux)
    result.append(v)
    return result
