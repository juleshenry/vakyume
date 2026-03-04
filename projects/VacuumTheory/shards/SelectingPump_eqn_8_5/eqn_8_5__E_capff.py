from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_5__Eff(
    self,
    actual_brake_horsepower: float,
    theoretical_adiabatic_horsepower: float,
    **kwargs,
):
    # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
    result = []
    Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
    result.append(Eff)
    return result
