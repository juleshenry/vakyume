from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_14__v_s(self, M: float, R: float, T: float, g_c: float, k: float, **kwargs):
    # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
    result = []
    v_s = sqrt(R * T * g_c * k / M)
    result.append(v_s)
    return result
