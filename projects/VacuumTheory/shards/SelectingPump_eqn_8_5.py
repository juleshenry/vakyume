from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SelectingPump:
    @kwasak_static
    def eqn_8_5(Eff=None, actual_brake_horsepower=None, theoretical_adiabatic_horsepower=None, **kwargs):
        return

    @staticmethod
    def eqn_8_5__Eff(actual_brake_horsepower: float, theoretical_adiabatic_horsepower: float):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        Eff = theoretical_adiabatic_horsepower/actual_brake_horsepower
        result.append(Eff)
        return result

    @staticmethod
    def eqn_8_5__actual_brake_horsepower(Eff: float, theoretical_adiabatic_horsepower: float):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        actual_brake_horsepower = theoretical_adiabatic_horsepower/Eff
        result.append(actual_brake_horsepower)
        return result

    @staticmethod
    def eqn_8_5__theoretical_adiabatic_horsepower(Eff: float, actual_brake_horsepower: float):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        theoretical_adiabatic_horsepower = Eff*actual_brake_horsepower
        result.append(theoretical_adiabatic_horsepower)
        return result

