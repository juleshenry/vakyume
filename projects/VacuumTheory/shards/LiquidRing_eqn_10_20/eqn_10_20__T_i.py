from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_20__T_i(self, P, S_0, S_p, T_e, p_0, p_c, p_s, **kwargs):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    result = []
    numerator = (S_0 / S_p) * ((P - p_0) * (460 + T_e)) * (P - p_c)
    denominator = P * (P - p_s) * (460 + T_e)
    term1 = pow(numerator / denominator, 1 / 3)
    term2 = ((P - p_0) * (460 + T_i) * (P - p_c)) / (P * (P - p_s) * (460 + T_e))
    result.append(
        (term1 ** (-5 / 3))
        * (
            (460 + T_i)
            * ((P - p_0) * (460 + T_e) * (P - p_c) / (P * (P - p_s) * (460 + T_e)))
        )
        ** (-2 / 3)
    )
    return [result]
