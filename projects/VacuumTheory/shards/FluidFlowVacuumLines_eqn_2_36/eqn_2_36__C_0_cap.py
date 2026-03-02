from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_36__C_0(self, C: float, F_t: float, **kwargs):
    # [.pyeqn] C = C_0 * F_t
    result = []
    C_0 = C/F_t
    result.append(C_0)
    return result
