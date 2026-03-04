from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_2__u_0(self, c: float, u_x: float, v: float, x: float, **kwargs):
    # [.pyeqn] u_x = u_0 * x / (1 - v * u_x / c ** 2)
    result = []
    u_0 = u_x * (c**2 - u_x * v) / (c**2 * x)
    result.append(u_0)
    return result
