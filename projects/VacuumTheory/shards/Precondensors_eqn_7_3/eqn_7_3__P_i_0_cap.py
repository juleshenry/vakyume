from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_3__P_i_0(epsilon_i: float, p_i: float, x_i: float, **kwargs):
    # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
    result = []
    P_i_0 = p_i/(epsilon_i*x_i)
    result.append(P_i_0)
    return result
