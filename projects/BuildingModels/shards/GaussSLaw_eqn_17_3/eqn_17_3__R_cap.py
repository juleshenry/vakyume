from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_17_3__R(self, E: float, Q: float, **kwargs):
    # [.pyeqn] E = Q / (4 * R ** 2)
    result = []
    R = -sqrt(Q / E) / 2
    result.append(R)
    R = sqrt(Q / E) / 2
    result.append(R)
    return result
