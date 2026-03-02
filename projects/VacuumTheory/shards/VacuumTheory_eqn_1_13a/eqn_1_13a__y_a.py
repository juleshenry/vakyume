from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_13a__y_a(self, n: float, n_a: float, **kwargs):
    # [.pyeqn] y_a = n_a / n
    result = []
    y_a = n_a/n
    result.append(y_a)
    return result
