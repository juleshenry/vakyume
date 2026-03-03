from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_8__P_2(self, P_1: float, adiabatic_power_watts: float, f: float, **kwargs):
    # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
    # P_2 = P_1 * ratio
    ratio = (12.0 * adiabatic_power_watts / f + 1.0) ** (1.0 / 0.286)
    P_2 = P_1 * ratio
    return [P_2]
