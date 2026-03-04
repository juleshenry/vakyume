from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_5_2b__x_2(self, K_1: float, K_2: float, x_1: float, y_1: float, y_2: float, **kwargs):
    # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
    result = []
    x_2 = K_1*x_1*y_2/(K_2*y_1)
    result.append(x_2)
    return result
