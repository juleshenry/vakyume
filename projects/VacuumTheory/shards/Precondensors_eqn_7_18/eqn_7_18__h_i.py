from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_7_18__h_i(
    self,
    D_0: float,
    D_LM: float,
    D_i: float,
    R_fi: float,
    R_fo: float,
    R_nc: float,
    U_0: float,
    h_c: float,
    k_w: float,
    x_w: float,
    **kwargs,
):
    # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
    result = []
    h_i = (
        -D_0
        * D_LM
        * U_0
        * h_c
        * k_w
        / (
            D_0 * D_LM * R_fi * U_0 * h_c * k_w
            + D_0 * D_i * U_0 * h_c * x_w
            + D_LM * D_i * R_fo * U_0 * h_c * k_w
            + D_LM * D_i * R_nc * U_0 * h_c * k_w
            + D_LM * D_i * U_0 * k_w
            - D_LM * D_i * h_c * k_w
        )
    )
    result.append(h_i)
    return result
