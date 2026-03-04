from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_6__P_1(
    self,
    M: float,
    P_2: float,
    R: float,
    T: float,
    adiabatic_hp: float,
    k: float,
    w: float,
    **kwargs,
):
    # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
    result = []
    P_1 = P_2 / (
        1980000 * M * adiabatic_hp / (R * T * w)
        - 1980000 * M * adiabatic_hp / (R * T * k * w)
        + 1
    ) ** (k / (k - 1))
    result.append(P_1)
    return result
