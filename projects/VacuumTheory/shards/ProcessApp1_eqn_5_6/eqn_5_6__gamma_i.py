from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_6__gamma_i(P_0_i: float, p_i: float, x_i: float, **kwargs):
    # [.pyeqn] p_i = x_i * gamma_i * P_0_i
    result = []
    gamma_i = p_i/(P_0_i*x_i)
    result.append(gamma_i)
    return result
