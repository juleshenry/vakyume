from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SelectingPump:
    @kwasak_static
    def eqn_8_8(P_1=None, P_2=None, adiabatic_power_watts=None, f=None, **kwargs):
        return

    @staticmethod
    def eqn_8_8__P_1(P_2: float, adiabatic_power_watts: float, f: float, **kwargs):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        # Error during Sympy solve: Sympy solve failed
        def func(P_1):
            # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
            return eval("(f / 12 * ((P_2 / x) ** 0.286 - 1)) - (adiabatic_power_watts)".replace('x', str(P_1)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_8_8__P_2(P_1: float, adiabatic_power_watts: float, f: float, **kwargs):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        # Error during Sympy solve: Sympy solve failed
        def func(P_2):
            # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
            return eval("(f / 12 * ((x / P_1) ** 0.286 - 1)) - (adiabatic_power_watts)".replace('x', str(P_2)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_8_8__adiabatic_power_watts(P_1: float, P_2: float, f: float, **kwargs):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        adiabatic_power_watts = 0.0833333333333333*f*((P_2/P_1)**(143/500) - 1.0)
        result.append(adiabatic_power_watts)
        return result

    @staticmethod
    def eqn_8_8__f(P_1: float, P_2: float, adiabatic_power_watts: float, **kwargs):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        f = 12.0*adiabatic_power_watts/((P_2/P_1)**0.286 - 1.0)
        result.append(f)
        return result

