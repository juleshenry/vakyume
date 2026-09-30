from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_7__epsilon_i(
    self,
    M: float,
    P: float,
    P_i_0: float,
    W_air: float,
    W_i: float,
    p_c: float,
    x_i: float,
    **kwargs,
):
    # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
    result = []
    epsilon_i = 29 * W_i * (P - p_c) / (M * P_i_0 * W_air * x_i)
    result.append(epsilon_i)
    return result
