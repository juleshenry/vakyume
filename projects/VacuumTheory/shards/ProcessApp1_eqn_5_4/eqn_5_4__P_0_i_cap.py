from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_4__P_0_i(self, P: float, x_i: float, y_i: float, **kwargs):
    # [.pyeqn] y_i * P = x_i * P_0_i
    result = []
    P_0_i = P*y_i/x_i
    result.append(P_0_i)
    return result
