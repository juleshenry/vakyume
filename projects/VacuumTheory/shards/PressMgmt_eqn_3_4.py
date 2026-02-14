from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class PressMgmt:
    @kwasak_static
    def eqn_3_4(KAPPA=None, P=None, V=None, **kwargs):
        return

    @staticmethod
    def eqn_3_4__KAPPA(P: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        KAPPA = P*V
        result.append(KAPPA)
        return result

    @staticmethod
    def eqn_3_4__P(KAPPA: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        P = KAPPA/V
        result.append(P)
        return result

    @staticmethod
    def eqn_3_4__V(KAPPA: float, P: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        V = KAPPA/P
        result.append(V)
        return result

