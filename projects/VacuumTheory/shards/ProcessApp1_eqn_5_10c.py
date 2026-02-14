from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_10c(D=None, L_0=None, R=None, **kwargs):
        return

    @staticmethod
    def eqn_5_10c__D(L_0: float, R: float, **kwargs):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        D = L_0/R
        result.append(D)
        return result

    @staticmethod
    def eqn_5_10c__L_0(D: float, R: float, **kwargs):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        L_0 = D*R
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10c__R(D: float, L_0: float, **kwargs):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        R = L_0/D
        result.append(R)
        return result

