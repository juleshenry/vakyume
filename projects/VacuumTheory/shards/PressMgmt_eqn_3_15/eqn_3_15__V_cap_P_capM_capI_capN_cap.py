from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_3_15__V_PMIN(self, **kwargs):
    # [.pyeqn] V_PMIN = 3.141592653589793 / 4
    result = []
    V_PMIN = 0.785398163397448
    result.append(V_PMIN)
    return result
