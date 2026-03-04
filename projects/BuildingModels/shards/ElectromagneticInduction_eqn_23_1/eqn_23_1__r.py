from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_23_1__r(self, B0: float, E: float, a: float, **kwargs):
    # [.pyeqn] E = B0 * a * r ** 2
    result = []
    r = -sqrt(E / (B0 * a))
    result.append(r)
    r = sqrt(E / (B0 * a))
    result.append(r)
    return result
