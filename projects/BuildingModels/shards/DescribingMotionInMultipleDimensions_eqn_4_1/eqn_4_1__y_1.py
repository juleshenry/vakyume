from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_4_1__y_1(
    self, t: float, v: float, x_1: float, x_2: float, y_2: float, **kwargs
):
    # [.pyeqn] v = (x_2 - x_1) / t + (y_2 - y_1) / t
    result = []
    y_1 = -t * v - x_1 + x_2 + y_2
    result.append(y_1)
    return result
