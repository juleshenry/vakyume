from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_4_1__x_2(
    self, t: float, v: float, x_1: float, y_1: float, y_2: float, **kwargs
):
    # [.pyeqn] v = (x_2 - x_1) / t + (y_2 - y_1) / t
    result = []
    x_2 = t * v + x_1 + y_1 - y_2
    result.append(x_2)
    return result
