from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_5_1__F(self, a: float, m: float, **kwargs):
    # [.pyeqn] F = m * a
    result = []
    F = a * m
    result.append(F)
    return result
