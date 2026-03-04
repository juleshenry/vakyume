from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_4_6__v_0(self, a: float, r: float, r_0: float, t: float, **kwargs):
    # [.pyeqn] r = r_0 + v_0 * t + 0.5 * a * t ** 2
    result = []
    v_0 = (-0.5 * a * t**2 + r - r_0) / t
    result.append(v_0)
    return result
