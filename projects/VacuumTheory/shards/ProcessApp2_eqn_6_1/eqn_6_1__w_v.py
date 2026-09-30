from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_6_1__w_v(
    self,
    T_1: float,
    T_2: float,
    T_R: float,
    c_p: float,
    del_h_v: float,
    w_1: float,
    w_2: float,
    **kwargs,
):
    # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
    result = []
    w_v = c_p * (T_1 * w_1 + T_2 * w_2 - T_R * w_1 - T_R * w_2) / del_h_v
    result.append(w_v)
    return result
