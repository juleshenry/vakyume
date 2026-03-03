from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_13b__y_a(self, P: float, p_a: float, **kwargs):
    # [.pyeqn] y_a = p_a / P
    result = []
    y_a = p_a/P
    result.append(y_a)
    return result
