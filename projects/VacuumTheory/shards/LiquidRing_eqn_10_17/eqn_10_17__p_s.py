from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_17__p_s(self, P: float, S_0: float, S_Th: float, p_0: float, **kwargs):
    # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
    result = []
    p_s = -P * (S_Th / S_0) ** (5 / 3) + P + p_0 * (S_Th / S_0) ** (5 / 3)
    result.append(p_s)
    p_s = (
        -P
        * (
            -0.5 * (S_Th / S_0) ** 0.333333333333333
            - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
        )
        ** 5
        + P
        + p_0
        * (
            -0.5 * (S_Th / S_0) ** 0.333333333333333
            - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
        )
        ** 5
    )
    result.append(p_s)
    p_s = (
        -P
        * (
            -0.5 * (S_Th / S_0) ** 0.333333333333333
            + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
        )
        ** 5
        + P
        + p_0
        * (
            -0.5 * (S_Th / S_0) ** 0.333333333333333
            + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
        )
        ** 5
    )
    result.append(p_s)
    return result
