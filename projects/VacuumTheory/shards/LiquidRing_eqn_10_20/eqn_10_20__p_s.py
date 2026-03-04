from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_20__p_s(self, P, S_0, S_p, T_e, T_i, p_0, p_c, **kwargs):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    result = []
    numerator = (S_0 / S_p) * ((P - p_0) * (460 + T_i) * (P - p_c))
    denominator = P * (P - (p_0 * p_s / (460 + T_i))) * (460 + T_e)
    p_s = ((numerator / denominator) ** (1 / 0.6) - (p_0 * p_c)) / (P - p_c)
    result.append(p_s)
    return [result]
