from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class PressMgmt:
    @kwasak_static
    def eqn_3_12(H_2=None, KAPPA_1=None, P=None, **kwargs):
        return

    @staticmethod
    def eqn_3_12__H_2(KAPPA_1: float, P: float):
        # [.pyeqn] P = KAPPA_1 * H_2 ** 2
        result = []
        H_2 = -sqrt(P/KAPPA_1)
        result.append(H_2)
        H_2 = sqrt(P/KAPPA_1)
        result.append(H_2)
        return result

    @staticmethod
    def eqn_3_12__KAPPA_1(H_2: float, P: float):
        # [.pyeqn] P = KAPPA_1 * H_2 ** 2
        result = []
        KAPPA_1 = P/H_2**2
        result.append(KAPPA_1)
        return result

    @staticmethod
    def eqn_3_12__P(H_2: float, KAPPA_1: float):
        # [.pyeqn] P = KAPPA_1 * H_2 ** 2
        result = []
        P = H_2**2*KAPPA_1
        result.append(P)
        return result

