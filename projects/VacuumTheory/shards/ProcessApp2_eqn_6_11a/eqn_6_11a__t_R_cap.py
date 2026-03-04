from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_6_11a__t_R(
    self,
    A_d: float,
    delta_T: float,
    delta_h_i: float,
    delta_m: float,
    h_d: float,
    m_b: float,
    **kwargs,
):
    # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
    result = []
    t_R = delta_h_i * delta_m * m_b / (A_d * delta_T * h_d)
    result.append(t_R)
    return result
