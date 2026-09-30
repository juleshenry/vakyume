from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_6_7__c_p(
    self,
    C_1: float,
    C_2: float,
    T_1: float,
    T_2: float,
    delta_h_c: float,
    delta_h_v: float,
    m_b: float,
    m_v: float,
    **kwargs,
):
    # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
    result = []
    c_p = (-C_1 * delta_h_c * m_b + C_2 * delta_h_c * m_b + delta_h_v * m_v) / (
        m_b * (T_1 - T_2)
    )
    result.append(c_p)
    return result
