from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_5_6__p_i(self, P_0_i: float, gamma_i: float, x_i: float, **kwargs):
    # [.pyeqn] p_i = x_i * gamma_i * P_0_i
    result = []
    p_i = P_0_i * gamma_i * x_i
    result.append(p_i)
    return result
