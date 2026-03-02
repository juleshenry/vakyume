from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_14__T(self, M: float, R: float, g_c: float, k: float, v_s: float, **kwargs):
    # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
    result = []
    T = M*v_s**2/(R*g_c*k)
    result.append(T)
    return result
