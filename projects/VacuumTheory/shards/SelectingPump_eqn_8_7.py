from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SelectingPump:
    @kwasak_static
    def eqn_8_7(P_1=None, P_2=None, adiabatic_hp=None, w=None, **kwargs):
        return

    @staticmethod
    def eqn_8_7__P_1(P_2: float, adiabatic_hp: float, w: float, **kwargs):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        # Error during Sympy solve: Sympy solve failed
        def func(P_1):
            # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
            return eval("((w / 20) * ((P_2 / x) ** 0.286 - 1)) - (adiabatic_hp)".replace('x', str(P_1)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_8_7__P_2(P_1: float, adiabatic_hp: float, w: float, **kwargs):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        # Error during Sympy solve: Sympy solve failed
        def func(P_2):
            # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
            return eval("((w / 20) * ((x / P_1) ** 0.286 - 1)) - (adiabatic_hp)".replace('x', str(P_2)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_8_7__adiabatic_hp(P_1: float, P_2: float, w: float, **kwargs):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        adiabatic_hp = 0.05*w*((P_2/P_1)**(143/500) - 1.0)
        result.append(adiabatic_hp)
        return result

    @staticmethod
    def eqn_8_7__w(P_1: float, P_2: float, adiabatic_hp: float, **kwargs):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        w = 20.0*adiabatic_hp/((P_2/P_1)**0.286 - 1.0)
        result.append(w)
        return result

