from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_13b__P(self, p_a: float, y_a: float, **kwargs):
    # [.pyeqn] y_a = p_a / P
    result = []
    P = p_a/y_a
    result.append(P)
    return result
