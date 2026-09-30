from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_6__R(
    self,
    M: float,
    P_1: float,
    P_2: float,
    T: float,
    adiabatic_hp: float,
    k: float,
    w: float,
    **kwargs,
):
    # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
    result = []
    R = (
        1980000
        * M
        * adiabatic_hp
        * (k - 1)
        / (T * k * w * ((P_2 / P_1) ** ((k - 1) / k) - 1))
    )
    result.append(R)
    return result
