from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_5_1__K_i(self, x_i: float, y_i: float, **kwargs):
    # [.pyeqn] K_i = y_i / x_i
    result = []
    K_i = y_i / x_i
    result.append(K_i)
    return result
