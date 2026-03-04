from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_5_1__a(self, F: float, m: float, **kwargs):
    # [.pyeqn] F = m * a
    result = []
    a = F / m
    result.append(a)
    return result
