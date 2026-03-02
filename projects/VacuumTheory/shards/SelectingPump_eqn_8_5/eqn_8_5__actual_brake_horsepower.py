from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_5__actual_brake_horsepower(self, Eff: float, theoretical_adiabatic_horsepower: float, **kwargs):
    # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
    result = []
    actual_brake_horsepower = theoretical_adiabatic_horsepower/Eff
    result.append(actual_brake_horsepower)
    return result
