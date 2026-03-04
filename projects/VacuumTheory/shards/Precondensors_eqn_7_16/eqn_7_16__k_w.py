from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_16__k_w(
    self,
    D_0: float,
    D_LM: float,
    D_i: float,
    R_f_0: float,
    R_fi: float,
    U_0: float,
    h_0: float,
    h_i: float,
    x_w: float,
    **kwargs,
):
    # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
    result = []
    k_w = (
        -D_0
        * D_i
        * U_0
        * h_0
        * h_i
        * x_w
        / (
            D_LM
            * (
                D_0 * R_fi * U_0 * h_0 * h_i
                + D_0 * U_0 * h_0
                + D_i * R_f_0 * U_0 * h_0 * h_i
                + D_i * U_0 * h_i
                - D_i * h_0 * h_i
            )
        )
    )
    result.append(k_w)
    return result
