from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_1__Abs_Pressure(BarometricPressure: float, Vacuum: float, **kwargs):
    # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
    result = []
    Abs_Pressure = BarometricPressure - Vacuum
    result.append(Abs_Pressure)
    return result
