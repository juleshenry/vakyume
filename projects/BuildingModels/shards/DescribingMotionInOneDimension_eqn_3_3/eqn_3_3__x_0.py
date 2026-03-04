from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_3__x_0(self, ax: float, t: float, v_0x: float, x: float, **kwargs):
    # [.pyeqn] x = x_0 + v_0x * t + 0.5 * ax * t ** 2
    result = []
    x_0 = -0.5 * ax * t**2 - t * v_0x + x
    result.append(x_0)
    return result
