from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_16__U_0(
    self,
    D_0: float,
    D_LM: float,
    D_i: float,
    R_f_0: float,
    R_fi: float,
    h_0: float,
    h_i: float,
    k_w: float,
    x_w: float,
    **kwargs,
):
    # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
    result = []
    U_0 = (
        D_LM
        * D_i
        * h_0
        * h_i
        * k_w
        / (
            D_0 * D_LM * R_fi * h_0 * h_i * k_w
            + D_0 * D_LM * h_0 * k_w
            + D_0 * D_i * h_0 * h_i * x_w
            + D_LM * D_i * R_f_0 * h_0 * h_i * k_w
            + D_LM * D_i * h_i * k_w
        )
    )
    result.append(U_0)
    return result
