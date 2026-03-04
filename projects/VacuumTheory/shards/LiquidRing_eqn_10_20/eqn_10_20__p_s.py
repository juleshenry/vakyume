from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_20__p_s(
    self,
    P: float,
    S_0: float,
    S_p: float,
    T_e: float,
    T_i: float,
    p_0: float,
    p_c: float,
    **kwargs,
):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    result = []
    p_s = (
        P**2 * T_e * (S_0 / S_p) ** (5 / 3)
        - P**2 * T_i
        + 460.0 * P**2 * (S_0 / S_p) ** (5 / 3)
        - 460.0 * P**2
        + P * T_i * p_0
        + P * T_i * p_c
        + 460.0 * P * p_0
        + 460.0 * P * p_c
        - T_i * p_0 * p_c
        - 460.0 * p_0 * p_c
    ) / (P * (S_0 / S_p) ** (5 / 3) * (T_e + 460.0))
    result.append(p_s)
    p_s = (
        P**2
        * T_e
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
        - P**2 * T_i
        + 460.0
        * P**2
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
        - 460.0 * P**2
        + P * T_i * p_0
        + P * T_i * p_c
        + 460.0 * P * p_0
        + 460.0 * P * p_c
        - T_i * p_0 * p_c
        - 460.0 * p_0 * p_c
    ) / (
        P
        * (T_e + 460.0)
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
    )
    result.append(p_s)
    p_s = (
        P**2
        * T_e
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
        - P**2 * T_i
        + 460.0
        * P**2
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
        - 460.0 * P**2
        + P * T_i * p_0
        + P * T_i * p_c
        + 460.0 * P * p_0
        + 460.0 * P * p_c
        - T_i * p_0 * p_c
        - 460.0 * p_0 * p_c
    ) / (
        P
        * (T_e + 460.0)
        * (
            -0.5 * (S_0 / S_p) ** 0.333333333333333
            + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
        )
        ** 5
    )
    result.append(p_s)
    return result
