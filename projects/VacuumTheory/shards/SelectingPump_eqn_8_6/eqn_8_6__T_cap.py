from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_6__T(
    self,
    M: float,
    P_1: float,
    P_2: float,
    R: float,
    adiabatic_hp: float,
    k: float,
    w: float,
    **kwargs,
):
    # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
    result = []
    T = (
        1980000
        * M
        * adiabatic_hp
        * (k - 1)
        / (R * k * w * ((P_2 / P_1) ** ((k - 1) / k) - 1))
    )
    result.append(T)
    return result
