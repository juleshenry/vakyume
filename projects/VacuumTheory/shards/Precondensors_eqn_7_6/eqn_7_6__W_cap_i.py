from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_6__W_i(
    self,
    M: float,
    P: float,
    P_i_0: float,
    W_air: float,
    p_c: float,
    x_i: float,
    **kwargs,
):
    # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
    result = []
    W_i = M * P_i_0 * W_air * x_i / (29 * (P - p_c))
    result.append(W_i)
    return result
