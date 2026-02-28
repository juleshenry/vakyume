from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_3__x_i(P_0_i: float, p_i: float, **kwargs):
    # [.pyeqn] p_i = x_i * P_0_i
    result = []
    x_i = p_i/P_0_i
    result.append(x_i)
    return result
