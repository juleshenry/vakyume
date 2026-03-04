from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_20__T_i(self, P, S_0, S_p, T_e, p_0, p_c, p_s, **kwargs):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    result = []
    T_i = (
        P**2 * T_e * (S_0 / S_p) ** (5 / 3)
        + 460.0 * P**2 * (S_0 / S_p) ** (5 / 3)
        - 460.0 * P**2
        - P * T_e * p_s * (S_0 / S_p) ** (5 / 3)
        + 460.0 * P * p_0
        + 460.0 * P * p_c
        - 460.0 * P * p_s * (S_0 / S_p) ** (5 / 3)
        - 460.0 * p_0 * p_c
    ) / (P**2 - P * p_0 - P * p_c + p_0 * p_c)
    result.append(T_i)
    T_i = (
        P**2
        * T_e
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
        + 460.0
        * P**2
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
        - 460.0 * P**2
        - P
        * T_e
        * p_s
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
        + 460.0 * P * p_0
        + 460.0 * P * p_c
        - 460.0
        * P
        * p_s
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
        - 460.0 * p_0 * p_c
    ) / (P**2 - P * p_0 - P * p_c + p_0 * p_c)
    result.append(T_i)
    T_i = (
        P**2
        * T_e
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
        + 460.0
        * P**2
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
        - 460.0 * P**2
        - P
        * T_e
        * p_s
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
        + 460.0 * P * p_0
        + 460.0 * P * p_c
        - 460.0
        * P
        * p_s
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
        - 460.0 * p_0 * p_c
    ) / (P**2 - P * p_0 - P * p_c + p_0 * p_c)
    result.append(T_i)
    return result
