from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_20__P(self, S_0, S_p, T_e, T_i, p_0, p_c, p_s, **kwargs):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    result = []
    numerator = (S_0 / S_p) * ((T_e + 460) * (p_s - P) * (P - p_c))
    denominator = ((P - p_0) * (460 + T_i)) * (P * (P - p_s) * (460 + T_e))
    P = pow(numerator / denominator, 1 / 0.6)
    result.append(P)
    return [result]
