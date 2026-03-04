from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_2__u_x(self, c: float, u_0: float, v: float, x: float, **kwargs):
    # [.pyeqn] u_x = u_0 * x / (1 - v * u_x / c ** 2)
    result = []
    u_x = c * (c - sqrt(c**2 - 4 * u_0 * v * x)) / (2 * v)
    result.append(u_x)
    u_x = c * (c + sqrt(c**2 - 4 * u_0 * v * x)) / (2 * v)
    result.append(u_x)
    return result
