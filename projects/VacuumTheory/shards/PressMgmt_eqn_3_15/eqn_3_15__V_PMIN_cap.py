from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_15__V_PMIN(self, **kwargs):
    # [.pyeqn] V_PMIN = 3.141592653589793 / 4
    result = []
    V_PMIN = 0.785398163397448
    result.append(V_PMIN)
    return result
