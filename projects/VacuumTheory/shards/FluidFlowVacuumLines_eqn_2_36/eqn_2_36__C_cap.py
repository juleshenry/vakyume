from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_36__C(self, C_0: float, F_t: float, **kwargs):
    # [.pyeqn] C = C_0 * F_t
    result = []
    C = C_0*F_t
    result.append(C)
    return result
