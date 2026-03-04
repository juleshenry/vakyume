from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_9_3__G(self, M: float, R: float, a: float, m: float, **kwargs):
    # [.pyeqn] a = G * M * m / R ** 2
    result = []
    G = R**2 * a / (M * m)
    result.append(G)
    return result
