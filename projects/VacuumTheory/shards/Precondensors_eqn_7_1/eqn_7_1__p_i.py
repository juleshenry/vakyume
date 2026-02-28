from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_1__p_i(P: float, y_i: float, **kwargs):
    # [.pyeqn] y_i = p_i / P
    result = []
    p_i = P*y_i
    result.append(p_i)
    return result
