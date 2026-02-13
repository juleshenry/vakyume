from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class PressMgmt:
    @kwasak_static
    def eqn_3_04(KAPPA: float = None, P: float = None, V: float = None, **kwargs):
        return

    @staticmethod
    def eqn_3_04__KAPPA(P: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        KAPPA = P * V
        result.append(KAPPA)
        return result

    @staticmethod
    def eqn_3_04__P(KAPPA: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        P = KAPPA / V
        result.append(P)
        return result

    @staticmethod
    def eqn_3_04__V(KAPPA: float, P: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        V = KAPPA / P
        result.append(V)
        return result


