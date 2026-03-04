from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_5_8__P_0_1(
    self, P_0_2: float, alpha_12: float, gamma_1: float, gamma_2: float, **kwargs
):
    # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
    result = []
    P_0_1 = P_0_2 * alpha_12 * gamma_2 / gamma_1
    result.append(P_0_1)
    return result
