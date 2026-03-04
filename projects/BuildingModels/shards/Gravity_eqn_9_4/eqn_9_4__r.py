from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_9_4__r(self, G: float, M: float, U: float, m: float, **kwargs):
    # [.pyeqn] U = - G * M * m / r
    result = []
    r = -G * M * m / U
    result.append(r)
    return result
