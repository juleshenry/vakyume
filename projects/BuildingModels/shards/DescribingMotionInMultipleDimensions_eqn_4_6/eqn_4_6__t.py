from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_4_6__t(self, a: float, r: float, r_0: float, v_0: float, **kwargs):
    # [.pyeqn] r = r_0 + v_0 * t + 0.5 * a * t ** 2
    result = []
    t = (-v_0 - sqrt(2.0 * a * r - 2.0 * a * r_0 + v_0**2)) / a
    result.append(t)
    t = (-v_0 + sqrt(2.0 * a * r - 2.0 * a * r_0 + v_0**2)) / a
    result.append(t)
    return result
