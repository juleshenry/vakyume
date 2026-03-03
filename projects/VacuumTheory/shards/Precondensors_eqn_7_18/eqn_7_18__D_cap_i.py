from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_7_18__D_i(
    self,
    D_0: float,
    D_LM: float,
    R_fi: float,
    R_fo: float,
    R_nc: float,
    U_0: float,
    h_c: float,
    h_i: float,
    k_w: float,
    x_w: float,
    **kwargs,
):
    # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
    result = []
    D_i = (
        -D_0
        * D_LM
        * U_0
        * h_c
        * k_w
        * (R_fi * h_i + 1)
        / (
            h_i
            * (
                D_0 * U_0 * h_c * x_w
                + D_LM * R_fo * U_0 * h_c * k_w
                + D_LM * R_nc * U_0 * h_c * k_w
                + D_LM * U_0 * k_w
                - D_LM * h_c * k_w
            )
        )
    )
    result.append(D_i)
    return result
