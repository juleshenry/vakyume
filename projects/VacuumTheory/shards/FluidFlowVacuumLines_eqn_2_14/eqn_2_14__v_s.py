from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_14__v_s(M: float, R: float, T: float, g_c: float, k: float, **kwargs):
    # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
    result = []
    v_s = sqrt(R*T*g_c*k/M)
    result.append(v_s)
    return result
