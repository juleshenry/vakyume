from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_1__BarometricPressure(self, Abs_Pressure: float, Vacuum: float, **kwargs):
    # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
    result = []
    BarometricPressure = Abs_Pressure + Vacuum
    result.append(BarometricPressure)
    return result
