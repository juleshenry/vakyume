from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_17__P(self, S_0: float, S_Th: float, p_0: float, p_s: float, **kwargs):
    # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
    result = []
    P = (p_0 * (S_Th / S_0) ** (5 / 3) - p_s) / ((S_Th / S_0) ** 1.66666666666667 - 1.0)
    result.append(P)
    P = (
        p_0
        * (
            -0.5 * (S_Th / S_0) ** 0.333333333333333
            - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
        )
        ** 5
        - p_s
    ) / (
        (
            -0.5 * (S_Th / S_0) ** 0.333333333333333
            - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
        )
        ** 5
        - 1.0
    )
    result.append(P)
    P = (
        p_0
        * (
            -0.5 * (S_Th / S_0) ** 0.333333333333333
            + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
        )
        ** 5
        - p_s
    ) / (
        (
            -0.5 * (S_Th / S_0) ** 0.333333333333333
            + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
        )
        ** 5
        - 1.0
    )
    result.append(P)
    return result
