from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_15__sum_R(self, U: float, **kwargs):
    # [.pyeqn] 1 / U = sum_R
    result = []
    sum_R = 1/U
    result.append(sum_R)
    return result
