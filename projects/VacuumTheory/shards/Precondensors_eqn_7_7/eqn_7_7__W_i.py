from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_7_7__W_i(
    self,
    M: float,
    P: float,
    P_i_0: float,
    W_air: float,
    epsilon_i: float,
    p_c: float,
    x_i: float,
    **kwargs,
):
    # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
    result = []
    W_i = M * P_i_0 * W_air * epsilon_i * x_i / (29 * (P - p_c))
    result.append(W_i)
    return result
