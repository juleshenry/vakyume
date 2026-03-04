from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_11_8__r(self, I: float, i: float, m: float, **kwargs):
    # [.pyeqn] I = (m * r ** 2) / i
    result = []
    r = -sqrt(I * i / m)
    result.append(r)
    r = sqrt(I * i / m)
    result.append(r)
    return result
