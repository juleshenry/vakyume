from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_20__p_c(self, P, S_0, S_p, T_e, T_i, p_0, p_s, **kwargs):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    result = []
    numerator = (S_0 / S_p) * ((P - p_0) * (460 + T_i) * (P - p_s))
    denominator = P * (P - T_e) * (460 + T_i)
    term1 = ((P - p_0) * (460 + T_i) * (P - p_s)) / (P * (P - p_s) * (460 + T_e))
    term2 = (numerator / denominator) ** (1 / 0.6)
    result.append((term1 - 1) * p_c)
    return [result[0]]
