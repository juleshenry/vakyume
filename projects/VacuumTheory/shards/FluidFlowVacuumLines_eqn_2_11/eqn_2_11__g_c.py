from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_11__g_c(self, D: float, L: float, f: float, h_r: float, v: float, **kwargs):
    # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
    result = []
    g_c = L*f*v**2/(2*D*h_r)
    result.append(g_c)
    return result
