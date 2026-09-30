from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_5_1__x_i(self, K_i: float, y_i: float, **kwargs):
    # [.pyeqn] K_i = y_i / x_i
    result = []
    x_i = y_i / K_i
    result.append(x_i)
    return result
