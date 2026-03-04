from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_19__P(
    self,
    S_Th: float,
    S_p: float,
    T_e: float,
    T_i: float,
    p_c: float,
    p_s: float,
    **kwargs,
):
    # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
    result = []
    P = (
        T_e * p_c * (S_p / S_Th) ** (5 / 3)
        - T_i * p_s
        + 460.0 * p_c * (S_p / S_Th) ** (5 / 3)
        - 460.0 * p_s
    ) / (
        T_e * (S_p / S_Th) ** 1.66666666666667
        - T_i
        + 460.0 * (S_p / S_Th) ** 1.66666666666667
        - 460.0
    )
    result.append(P)
    P = (
        T_e
        * p_c
        * (
            -0.5 * (S_p / S_Th) ** 0.333333333333333
            - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
        )
        ** 5
        - T_i * p_s
        + 460.0
        * p_c
        * (
            -0.5 * (S_p / S_Th) ** 0.333333333333333
            - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
        )
        ** 5
        - 460.0 * p_s
    ) / (
        T_e
        * (
            -0.5 * (S_p / S_Th) ** 0.333333333333333
            - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
        )
        ** 5
        - T_i
        + 460.0
        * (
            -0.5 * (S_p / S_Th) ** 0.333333333333333
            - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
        )
        ** 5
        - 460.0
    )
    result.append(P)
    P = (
        T_e
        * p_c
        * (
            -0.5 * (S_p / S_Th) ** 0.333333333333333
            + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
        )
        ** 5
        - T_i * p_s
        + 460.0
        * p_c
        * (
            -0.5 * (S_p / S_Th) ** 0.333333333333333
            + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
        )
        ** 5
        - 460.0 * p_s
    ) / (
        T_e
        * (
            -0.5 * (S_p / S_Th) ** 0.333333333333333
            + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
        )
        ** 5
        - T_i
        + 460.0
        * (
            -0.5 * (S_p / S_Th) ** 0.333333333333333
            + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
        )
        ** 5
        - 460.0
    )
    result.append(P)
    return result
