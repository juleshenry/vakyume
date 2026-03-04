from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_8__P_2(self, P_1, adiabatic_power_watts, f, **kwargs):
    # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
    result = []
    numerator = (adiabatic_power_watts * 12) + 1
    denominator = pow(P_1, 0.286)
    P_2 = numerator / denominator
    result.append(P_2)
    return [result]
