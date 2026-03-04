from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_9_4__m(self, G: float, M: float, U: float, r: float, **kwargs):
    # [.pyeqn] U = - G * M * m / r
    result = []
    m = -U * r / (G * M)
    result.append(m)
    return result
