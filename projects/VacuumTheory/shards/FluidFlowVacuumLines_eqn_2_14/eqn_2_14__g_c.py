from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_14__g_c(self, M: float, R: float, T: float, k: float, v_s: float, **kwargs):
    # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
    result = []
    g_c = M*v_s**2/(R*T*k)
    result.append(g_c)
    return result
