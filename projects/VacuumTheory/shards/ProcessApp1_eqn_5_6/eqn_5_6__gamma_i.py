from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_5_6__gamma_i(self, P_0_i: float, p_i: float, x_i: float, **kwargs):
    # [.pyeqn] p_i = x_i * gamma_i * P_0_i
    result = []
    gamma_i = p_i / (P_0_i * x_i)
    result.append(gamma_i)
    return result
