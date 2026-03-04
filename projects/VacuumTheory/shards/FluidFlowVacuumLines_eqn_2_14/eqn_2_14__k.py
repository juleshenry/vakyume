from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_14__k(self, M: float, R: float, T: float, g_c: float, v_s: float, **kwargs):
    # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
    result = []
    k = M * v_s**2 / (R * T * g_c)
    result.append(k)
    return result
