from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_6__M(
    self,
    P_1: float,
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
    M = (
        R
        * T
        * k
        * w
        * ((P_2 / P_1) ** ((k - 1) / k) - 1)
        / (1980000 * adiabatic_hp * (k - 1))
    )
    result.append(M)
    return result
