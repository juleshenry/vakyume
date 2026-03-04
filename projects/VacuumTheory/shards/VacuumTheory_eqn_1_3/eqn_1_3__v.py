from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_1_3__v(self, T: float, k: float, m: float, **kwargs):
    # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
    result = []
    v = -1.73205080756888 * sqrt(T * k / m)
    result.append(v)
    v = 1.73205080756888 * sqrt(T * k / m)
    result.append(v)
    return result
