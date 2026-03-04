from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_4_6__r(self, a: float, r_0: float, t: float, v_0: float, **kwargs):
    # [.pyeqn] r = r_0 + v_0 * t + 0.5 * a * t ** 2
    result = []
    r = 0.5 * a * t**2 + r_0 + t * v_0
    result.append(r)
    return result
