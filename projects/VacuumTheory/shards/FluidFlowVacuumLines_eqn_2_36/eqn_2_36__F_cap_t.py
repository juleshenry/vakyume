from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_36__F_t(self, C: float, C_0: float, **kwargs):
    # [.pyeqn] C = C_0 * F_t
    result = []
    F_t = C / C_0
    result.append(F_t)
    return result
