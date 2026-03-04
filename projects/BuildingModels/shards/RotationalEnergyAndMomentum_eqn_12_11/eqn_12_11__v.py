from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_12_11__v(self, I: float, L: float, **kwargs):
    # [.pyeqn] L = (1 / 2) * I * v ** 2
    result = []
    v = -sqrt(2) * sqrt(L / I)
    result.append(v)
    v = sqrt(2) * sqrt(L / I)
    result.append(v)
    return result
