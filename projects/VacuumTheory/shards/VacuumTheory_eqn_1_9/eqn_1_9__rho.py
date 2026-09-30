from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_1_9__rho(self, M: float, P: float, R: float, T: float, **kwargs):
    # [.pyeqn] rho = P * M / (R * T)
    result = []
    rho = M * P / (R * T)
    result.append(rho)
    return result
