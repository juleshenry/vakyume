from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_11__D(self, L: float, f: float, g_c: float, h_r: float, v: float, **kwargs):
    # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
    result = []
    D = L*f*v**2/(2*g_c*h_r)
    result.append(D)
    return result
