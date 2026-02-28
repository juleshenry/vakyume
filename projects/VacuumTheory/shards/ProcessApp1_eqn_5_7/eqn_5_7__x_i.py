from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_7__x_i(P: float, P_0_i: float, gamma_i: float, y_i: float, **kwargs):
    # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
    result = []
    x_i = P*y_i/(P_0_i*gamma_i)
    result.append(x_i)
    return result
