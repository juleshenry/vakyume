from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_1__K_i(x_i: float, y_i: float, **kwargs):
    # [.pyeqn] K_i = y_i / x_i
    result = []
    K_i = y_i/x_i
    result.append(K_i)
    return result
