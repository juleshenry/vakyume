from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_2__x_i(self, P_i_0: float, p_i: float, **kwargs):
    # [.pyeqn] p_i = x_i * P_i_0
    result = []
    x_i = p_i/P_i_0
    result.append(x_i)
    return result
