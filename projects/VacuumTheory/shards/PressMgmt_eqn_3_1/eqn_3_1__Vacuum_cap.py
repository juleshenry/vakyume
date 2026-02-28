from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_1__Vacuum(Abs_Pressure: float, BarometricPressure: float, **kwargs):
    # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
    result = []
    Vacuum = -Abs_Pressure + BarometricPressure
    result.append(Vacuum)
    return result
