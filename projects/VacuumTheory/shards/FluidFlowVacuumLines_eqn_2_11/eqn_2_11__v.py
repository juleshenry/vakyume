from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_11__v(self, D: float, L: float, f: float, g_c: float, h_r: float, **kwargs):
    # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
    result = []
    v = -sqrt(2)*sqrt(D*g_c*h_r/(L*f))
    result.append(v)
    v = sqrt(2)*sqrt(D*g_c*h_r/(L*f))
    result.append(v)
    return result
