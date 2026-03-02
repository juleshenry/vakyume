from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_7__y_i(self, P: float, P_0_i: float, gamma_i: float, x_i: float, **kwargs):
    # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
    result = []
    y_i = P_0_i*gamma_i*x_i/P
    result.append(y_i)
    return result
