from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_10b(L_0=None, R=None, V_1=None, **kwargs):
        return

    @staticmethod
    def eqn_5_10b__L_0(R: float, V_1: float, **kwargs):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        L_0 = R*V_1/(R + 1)
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10b__R(L_0: float, V_1: float, **kwargs):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        R = -L_0/(L_0 - V_1)
        result.append(R)
        return result

    @staticmethod
    def eqn_5_10b__V_1(L_0: float, R: float, **kwargs):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        V_1 = L_0 + L_0/R
        result.append(V_1)
        return result

