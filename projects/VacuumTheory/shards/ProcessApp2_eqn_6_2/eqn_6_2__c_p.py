from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_2__c_p(Q_v: float, T_1: float, T_2: float, T_R: float, w_1: float, w_2: float, **kwargs):
    # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
    result = []
    c_p = 12000*Q_v/(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)
    result.append(c_p)
    return result
