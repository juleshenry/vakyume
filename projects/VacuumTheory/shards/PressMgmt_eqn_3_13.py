from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class PressMgmt:
    @kwasak_static
    def eqn_3_13(H_1=None, H_2=None, KAPPA_2=None, P=None, **kwargs):
        return

    @staticmethod
    def eqn_3_13__H_1(H_2: float, KAPPA_2: float, P: float, **kwargs):
        # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
        result = []
        H_1 = H_2 - P/KAPPA_2
        result.append(H_1)
        return result

    @staticmethod
    def eqn_3_13__H_2(H_1: float, KAPPA_2: float, P: float, **kwargs):
        # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
        result = []
        H_2 = H_1 + P/KAPPA_2
        result.append(H_2)
        return result

    @staticmethod
    def eqn_3_13__KAPPA_2(H_1: float, H_2: float, P: float, **kwargs):
        # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
        result = []
        KAPPA_2 = -P/(H_1 - H_2)
        result.append(KAPPA_2)
        return result

    @staticmethod
    def eqn_3_13__P(H_1: float, H_2: float, KAPPA_2: float, **kwargs):
        # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
        result = []
        P = KAPPA_2*(-H_1 + H_2)
        result.append(P)
        return result

