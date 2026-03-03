from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_11__h_r(self, D: float, L: float, f: float, g_c: float, v: float, **kwargs):
    # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
    result = []
    h_r = L*f*v**2/(2*D*g_c)
    result.append(h_r)
    return result
