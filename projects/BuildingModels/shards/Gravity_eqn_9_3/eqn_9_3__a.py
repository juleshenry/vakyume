from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_9_3__a(self, G: float, M: float, R: float, m: float, **kwargs):
    # [.pyeqn] a = G * M * m / R ** 2
    result = []
    a = G * M * m / R**2
    result.append(a)
    return result
