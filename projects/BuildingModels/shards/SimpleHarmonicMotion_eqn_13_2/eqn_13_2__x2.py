from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_13_2__x2(self, x0: float, x1: float, **kwargs):
    # [.pyeqn] 1 * x1 + 2 * x2 = (1 + 2) * x0
    result = []
    x2 = 3 * x0 / 2 - x1 / 2
    result.append(x2)
    return result
