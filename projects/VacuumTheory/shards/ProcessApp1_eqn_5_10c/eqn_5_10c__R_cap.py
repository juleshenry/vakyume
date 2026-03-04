from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_5_10c__R(self, D: float, L_0: float, **kwargs):
    # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
    result = []
    R = L_0 / D
    result.append(R)
    return result
