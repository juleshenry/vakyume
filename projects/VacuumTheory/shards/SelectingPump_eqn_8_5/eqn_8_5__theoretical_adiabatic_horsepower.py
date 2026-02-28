from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_5__theoretical_adiabatic_horsepower(Eff: float, actual_brake_horsepower: float, **kwargs):
    # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
    result = []
    theoretical_adiabatic_horsepower = Eff*actual_brake_horsepower
    result.append(theoretical_adiabatic_horsepower)
    return result
