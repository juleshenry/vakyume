from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_8_5__Eff_cap import eqn_8_5__Eff
from .eqn_8_5__actual_brake_horsepower import eqn_8_5__actual_brake_horsepower
from .eqn_8_5__theoretical_adiabatic_horsepower import eqn_8_5__theoretical_adiabatic_horsepower

class SelectingPump:
    eqn_8_5__Eff = staticmethod(eqn_8_5__Eff)
    eqn_8_5__actual_brake_horsepower = staticmethod(eqn_8_5__actual_brake_horsepower)
    eqn_8_5__theoretical_adiabatic_horsepower = staticmethod(eqn_8_5__theoretical_adiabatic_horsepower)

    @kwasak_static
    def eqn_8_5(Eff=None, actual_brake_horsepower=None, theoretical_adiabatic_horsepower=None, **kwargs):
        return
