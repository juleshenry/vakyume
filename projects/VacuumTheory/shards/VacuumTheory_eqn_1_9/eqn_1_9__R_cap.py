from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_1_9__R(self, M: float, P: float, T: float, rho: float, **kwargs):
    # [.pyeqn] rho = P * M / (R * T)
    result = []
    R = M * P / (T * rho)
    result.append(R)
    return result
