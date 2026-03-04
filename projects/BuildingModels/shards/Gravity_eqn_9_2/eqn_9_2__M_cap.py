from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_9_2__M(self, G: float, R: float, T: float, pi: float, **kwargs):
    # [.pyeqn] T ** 2 = G * M * R ** 3 / (4 * pi ** 2)
    result = []
    M = 4 * T**2 * pi**2 / (G * R**3)
    result.append(M)
    return result
