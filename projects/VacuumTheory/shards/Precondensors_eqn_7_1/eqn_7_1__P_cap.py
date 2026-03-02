from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_1__P(self, p_i: float, y_i: float, **kwargs):
    # [.pyeqn] y_i = p_i / P
    result = []
    P = p_i/y_i
    result.append(P)
    return result
