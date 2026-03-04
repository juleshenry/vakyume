from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_6_1__F_N(self, F: float, **kwargs):
    # [.pyeqn] F = F_N
    result = []
    F_N = F
    result.append(F_N)
    return result
