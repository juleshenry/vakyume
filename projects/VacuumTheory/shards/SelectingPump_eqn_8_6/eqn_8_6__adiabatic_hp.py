from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_6__adiabatic_hp(self, M: float, P_1: float, P_2: float, R: float, T: float, k: float, w: float, **kwargs):
    # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
    result = []
    adiabatic_hp = R*T*k*w*((P_2/P_1)**((k - 1)/k) - 1)/(1980000*M*(k - 1))
    result.append(adiabatic_hp)
    return result
