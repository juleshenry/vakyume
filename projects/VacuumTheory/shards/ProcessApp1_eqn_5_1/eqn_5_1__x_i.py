from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_1__x_i(self, K_i: float, y_i: float, **kwargs):
    # [.pyeqn] K_i = y_i / x_i
    result = []
    x_i = y_i/K_i
    result.append(x_i)
    return result
