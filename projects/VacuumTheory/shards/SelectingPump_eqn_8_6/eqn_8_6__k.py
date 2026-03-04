from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_6__k(self, M, P_1, P_2, R, T, adiabatic_hp, w, **kwargs):
    # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
    result = []
    k = (
        M
        * 1980000
        * adiabatic_hp
        * (k - 1)
        / (R * T * w * ((P_2 / P_1) ** ((k - 1) / k) - 1))
    )
    if k == 0:
        return [None]
    result.append(k)
    return [result]
