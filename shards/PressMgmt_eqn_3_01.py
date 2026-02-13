from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class PressMgmt:
    @kwasak_static
    def eqn_3_01(
        Abs_Pressure: float = None,
        BarometricPressure: float = None,
        Vacuum: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_3_01__Abs_Pressure(BarometricPressure: float, Vacuum: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Abs_Pressure = BarometricPressure - Vacuum
        result.append(Abs_Pressure)
        return result

    @staticmethod
    def eqn_3_01__BarometricPressure(Abs_Pressure: float, Vacuum: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        BarometricPressure = Abs_Pressure + Vacuum
        result.append(BarometricPressure)
        return result

    @staticmethod
    def eqn_3_01__Vacuum(Abs_Pressure: float, BarometricPressure: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Vacuum = -Abs_Pressure + BarometricPressure
        result.append(Vacuum)
        return result


