from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_9_2__pi(self, G: float, M: float, R: float, T: float, **kwargs):
    # [.pyeqn] T ** 2 = G * M * R ** 3 / (4 * pi ** 2)
    result = []
    pi = -sqrt(G * M * R**3) / (2 * T)
    result.append(pi)
    pi = sqrt(G * M * R**3) / (2 * T)
    result.append(pi)
    return result
