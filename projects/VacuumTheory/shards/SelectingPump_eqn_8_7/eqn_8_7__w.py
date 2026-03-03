from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_7__w(self, P_1: float, P_2: float, adiabatic_hp: float, **kwargs):
    # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
    result = []
    w = 20.0*adiabatic_hp/((P_2/P_1)**0.286 - 1.0)
    result.append(w)
    return result
