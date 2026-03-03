from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_2__Q_v(self, T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float, **kwargs):
    # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
    result = []
    Q_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/12000
    result.append(Q_v)
    return result
