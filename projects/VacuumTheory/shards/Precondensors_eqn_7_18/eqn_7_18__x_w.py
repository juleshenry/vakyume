from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_18__x_w(
    self,
    D_0: float,
    D_LM: float,
    D_i: float,
    R_fi: float,
    R_fo: float,
    R_nc: float,
    U_0: float,
    h_c: float,
    h_i: float,
    k_w: float,
    **kwargs,
):
    # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
    result = []
    x_w = (
        -D_LM * R_fi * k_w / D_i
        - D_LM * k_w / (D_i * h_i)
        - D_LM * R_fo * k_w / D_0
        - D_LM * R_nc * k_w / D_0
        - D_LM * k_w / (D_0 * h_c)
        + D_LM * k_w / (D_0 * U_0)
    )
    result.append(x_w)
    return result
