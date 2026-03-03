from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_8__adiabatic_power_watts(self, P_1: float, P_2: float, f: float, **kwargs):
    # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
    result = []
    adiabatic_power_watts = 0.0833333333333333 * f * ((P_2 / P_1) ** (143 / 500) - 1.0)
    result.append(adiabatic_power_watts)
    return result
