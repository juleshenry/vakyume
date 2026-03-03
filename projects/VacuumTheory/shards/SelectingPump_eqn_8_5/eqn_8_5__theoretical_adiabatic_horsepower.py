from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_5__theoretical_adiabatic_horsepower(
    self, Eff: float, actual_brake_horsepower: float, **kwargs
):
    # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
    result = []
    theoretical_adiabatic_horsepower = Eff * actual_brake_horsepower
    result.append(theoretical_adiabatic_horsepower)
    return result
