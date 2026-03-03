from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_2__T_2(self, Q_v: float, T_1: float, T_R: float, c_p: float, w_1: float, w_2: float, **kwargs):
    # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
    result = []
    T_2 = (12000*Q_v + T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R))/(c_p*w_2)
    result.append(T_2)
    return result
