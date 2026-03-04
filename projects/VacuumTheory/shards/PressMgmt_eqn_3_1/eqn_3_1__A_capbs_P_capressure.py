from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_3_1__Abs_Pressure(self, BarometricPressure: float, Vacuum: float, **kwargs):
    # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
    result = []
    Abs_Pressure = BarometricPressure - Vacuum
    result.append(Abs_Pressure)
    return result
