from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_6_1__T_2(self, T_1: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float, **kwargs):
    # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
    result = []
    T_2 = (T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R) + del_h_v*w_v)/(c_p*w_2)
    result.append(T_2)
    return result
