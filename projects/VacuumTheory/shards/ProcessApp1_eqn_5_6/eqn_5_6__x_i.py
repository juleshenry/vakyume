from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_5_6__x_i(self, P_0_i: float, gamma_i: float, p_i: float, **kwargs):
    # [.pyeqn] p_i = x_i * gamma_i * P_0_i
    result = []
    x_i = p_i / (P_0_i * gamma_i)
    result.append(x_i)
    return result
