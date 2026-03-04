from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_6_8__T_2(
    self,
    C_1: float,
    C_2: float,
    T_1: float,
    c_p: float,
    delta_h_c: float,
    delta_h_v: float,
    delta_t: float,
    m_b: float,
    w_v: float,
    **kwargs,
):
    # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
    result = []
    T_2 = (
        T_1 * c_p * m_b + delta_h_c * m_b * (C_1 - C_2) - delta_h_v * delta_t * w_v
    ) / (c_p * m_b)
    result.append(T_2)
    return result
