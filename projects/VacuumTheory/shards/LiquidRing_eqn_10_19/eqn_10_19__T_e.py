from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_19__T_e(
    self,
    P: float,
    S_Th: float,
    S_p: float,
    T_i: float,
    p_c: float,
    p_s: float,
    **kwargs,
):
    # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
    result = []
    T_e = (
        P * T_i
        - 460.0 * P * (S_p / S_Th) ** (5 / 3)
        + 460.0 * P
        - T_i * p_s
        + 460.0 * p_c * (S_p / S_Th) ** (5 / 3)
        - 460.0 * p_s
    ) / ((S_p / S_Th) ** (5 / 3) * (P - p_c))
    result.append(T_e)
    T_e = (
        P * T_i
        - 460.0
        * P
        * (
            -0.5 * (S_p / S_Th) ** 0.333333333333333
            - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
        )
        ** 5
        + 460.0 * P
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
        (P - p_c)
        * (
            -0.5 * (S_p / S_Th) ** 0.333333333333333
            - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
        )
        ** 5
    )
    result.append(T_e)
    T_e = (
        P * T_i
        - 460.0
        * P
        * (
            -0.5 * (S_p / S_Th) ** 0.333333333333333
            + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
        )
        ** 5
        + 460.0 * P
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
        (P - p_c)
        * (
            -0.5 * (S_p / S_Th) ** 0.333333333333333
            + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
        )
        ** 5
    )
    result.append(T_e)
    return result
