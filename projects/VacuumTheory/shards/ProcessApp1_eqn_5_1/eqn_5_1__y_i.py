from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_5_1__y_i(self, K_i: float, x_i: float, **kwargs):
    # [.pyeqn] K_i = y_i / x_i
    result = []
    y_i = K_i * x_i
    result.append(y_i)
    return result
