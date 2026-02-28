from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_3_1__Abs_Pressure_cap import eqn_3_1__Abs_Pressure
from .eqn_3_1__BarometricPressure_cap import eqn_3_1__BarometricPressure
from .eqn_3_1__Vacuum_cap import eqn_3_1__Vacuum

class PressMgmt:
    eqn_3_1__Abs_Pressure = staticmethod(eqn_3_1__Abs_Pressure)
    eqn_3_1__BarometricPressure = staticmethod(eqn_3_1__BarometricPressure)
    eqn_3_1__Vacuum = staticmethod(eqn_3_1__Vacuum)

    @kwasak_static
    def eqn_3_1(Abs_Pressure=None, BarometricPressure=None, Vacuum=None, **kwargs):
        return
