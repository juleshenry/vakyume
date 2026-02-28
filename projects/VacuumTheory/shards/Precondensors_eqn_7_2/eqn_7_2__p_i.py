from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_2__p_i(P_i_0: float, x_i: float, **kwargs):
    # [.pyeqn] p_i = x_i * P_i_0
    result = []
    p_i = P_i_0*x_i
    result.append(p_i)
    return result
