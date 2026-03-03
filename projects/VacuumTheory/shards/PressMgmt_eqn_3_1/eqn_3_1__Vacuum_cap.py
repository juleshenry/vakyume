from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_1__Vacuum(self, Abs_Pressure: float, BarometricPressure: float, **kwargs):
    # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
    result = []
    Vacuum = -Abs_Pressure + BarometricPressure
    result.append(Vacuum)
    return result
