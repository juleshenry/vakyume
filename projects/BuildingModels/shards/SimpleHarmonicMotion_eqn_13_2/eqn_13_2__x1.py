from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_13_2__x1(self, x0: float, x2: float, **kwargs):
    # [.pyeqn] 1 * x1 + 2 * x2 = (1 + 2) * x0
    result = []
    x1 = 3 * x0 - 2 * x2
    result.append(x1)
    return result
