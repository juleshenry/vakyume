from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_18__R_fo(
    self,
    D_0: float,
    D_LM: float,
    D_i: float,
    R_fi: float,
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
    R_fo = (
        -D_0 * R_fi / D_i
        - D_0 / (D_i * h_i)
        - D_0 * x_w / (D_LM * k_w)
        - R_nc
        - 1 / h_c
        + 1 / U_0
    )
    result.append(R_fo)
    return result
