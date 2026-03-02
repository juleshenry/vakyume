from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_36__F_t(self, C: float, C_0: float, **kwargs):
    # [.pyeqn] C = C_0 * F_t
    result = []
    F_t = C/C_0
    result.append(F_t)
    return result
