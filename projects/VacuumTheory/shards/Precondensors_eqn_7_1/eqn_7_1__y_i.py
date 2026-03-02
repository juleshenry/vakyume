from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_1__y_i(self, P: float, p_i: float, **kwargs):
    # [.pyeqn] y_i = p_i / P
    result = []
    y_i = p_i/P
    result.append(y_i)
    return result
