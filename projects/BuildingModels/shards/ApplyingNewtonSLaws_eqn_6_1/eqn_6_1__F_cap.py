from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_6_1__F(self, F_N: float, **kwargs):
    # [.pyeqn] F = F_N
    result = []
    F = F_N
    result.append(F)
    return result
