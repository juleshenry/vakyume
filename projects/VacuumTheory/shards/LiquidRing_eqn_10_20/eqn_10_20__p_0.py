from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_20__p_0(self, P, S_0, S_p, T_e, T_i, p_c, p_s, **kwargs):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    result = []
    p_0 = (
        P
        * (
            -P * T_e * (S_0 / S_p) ** (5 / 3)
            + P * T_i
            - 460.0 * P * (S_0 / S_p) ** (5 / 3)
            + 460.0 * P
            + T_e * p_s * (S_0 / S_p) ** (5 / 3)
            - T_i * p_c
            - 460.0 * p_c
            + 460.0 * p_s * (S_0 / S_p) ** (5 / 3)
        )
        / (P * T_i + 460.0 * P - T_i * p_c - 460.0 * p_c)
    )
    result.append(p_0)
    p_0 = (
        P
        * (
            -P
            * T_e
            * (
                -0.5 * (S_0 / S_p) ** 0.333333333333333
                - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
            )
            ** 5
            + P * T_i
            - 460.0
            * P
            * (
                -0.5 * (S_0 / S_p) ** 0.333333333333333
                - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
            )
            ** 5
            + 460.0 * P
            + T_e
            * p_s
            * (
                -0.5 * (S_0 / S_p) ** 0.333333333333333
                - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
            )
            ** 5
            - T_i * p_c
            - 460.0 * p_c
            + 460.0
            * p_s
            * (
                -0.5 * (S_0 / S_p) ** 0.333333333333333
                - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
            )
            ** 5
        )
        / (P * T_i + 460.0 * P - T_i * p_c - 460.0 * p_c)
    )
    result.append(p_0)
    p_0 = (
        P
        * (
            -P
            * T_e
            * (
                -0.5 * (S_0 / S_p) ** 0.333333333333333
                + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
            )
            ** 5
            + P * T_i
            - 460.0
            * P
            * (
                -0.5 * (S_0 / S_p) ** 0.333333333333333
                + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
            )
            ** 5
            + 460.0 * P
            + T_e
            * p_s
            * (
                -0.5 * (S_0 / S_p) ** 0.333333333333333
                + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
            )
            ** 5
            - T_i * p_c
            - 460.0 * p_c
            + 460.0
            * p_s
            * (
                -0.5 * (S_0 / S_p) ** 0.333333333333333
                + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
            )
            ** 5
        )
        / (P * T_i + 460.0 * P - T_i * p_c - 460.0 * p_c)
    )
    result.append(p_0)
    return result
