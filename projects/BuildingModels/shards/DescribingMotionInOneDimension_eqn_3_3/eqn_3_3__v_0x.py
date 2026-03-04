from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_3__v_0x(self, ax: float, t: float, x: float, x_0: float, **kwargs):
    # [.pyeqn] x = x_0 + v_0x * t + 0.5 * ax * t ** 2
    result = []
    v_0x = (-0.5 * ax * t**2 + x - x_0) / t
    result.append(v_0x)
    return result
