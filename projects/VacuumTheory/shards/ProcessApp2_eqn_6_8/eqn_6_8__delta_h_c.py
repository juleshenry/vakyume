from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_6_8__delta_h_c(
    self,
    C_1: float,
    C_2: float,
    T_1: float,
    T_2: float,
    c_p: float,
    delta_h_v: float,
    delta_t: float,
    m_b: float,
    w_v: float,
    **kwargs,
):
    # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
    result = []
    delta_h_c = (-T_1 * c_p * m_b + T_2 * c_p * m_b + delta_h_v * delta_t * w_v) / (
        m_b * (C_1 - C_2)
    )
    result.append(delta_h_c)
    return result
