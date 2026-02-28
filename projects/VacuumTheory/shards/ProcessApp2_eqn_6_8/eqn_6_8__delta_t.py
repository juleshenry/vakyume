from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_8__delta_t(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, w_v: float, **kwargs):
    # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
    result = []
    delta_t = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*w_v)
    result.append(delta_t)
    return result
