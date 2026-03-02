from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_6__w(self, M: float, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float, **kwargs):
    # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
    result = []
    w = 1980000*M*adiabatic_hp*(k - 1)/(R*T*k*((P_2/P_1)**((k - 1)/k) - 1))
    result.append(w)
    return result
