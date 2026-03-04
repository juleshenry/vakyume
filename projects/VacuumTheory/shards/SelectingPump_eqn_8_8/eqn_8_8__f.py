from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_8_8__f(self, P_1: float, P_2: float, adiabatic_power_watts: float, **kwargs):
    # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
    result = []
    f = 12.0*adiabatic_power_watts/((P_2/P_1)**0.286 - 1.0)
    result.append(f)
    return result
