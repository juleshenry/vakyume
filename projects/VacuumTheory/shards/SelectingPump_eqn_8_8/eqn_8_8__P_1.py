from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_8__P_1(self, P_2, adiabatic_power_watts, f, **kwargs):
    # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
    return (f / (adiabatic_power_watts + 1.0)) ** (500 / 143) * P_2 ** (-1 / 0.286)
