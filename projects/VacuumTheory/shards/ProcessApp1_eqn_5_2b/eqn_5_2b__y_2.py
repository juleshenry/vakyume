from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_2b__y_2(K_1: float, K_2: float, x_1: float, x_2: float, y_1: float, **kwargs):
    # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
    result = []
    y_2 = K_2*x_2*y_1/(K_1*x_1)
    result.append(y_2)
    return result
