from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class PressMgmt:
    @kwasak_static
    def eqn_3_3(H_1=None, H_2=None, P=None, P_P=None, **kwargs):
        return

    @staticmethod
    def eqn_3_3__H_1(H_2: float, P: float, P_P: float, **kwargs):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_1 = H_2 + P - P_P
        result.append(H_1)
        return result

    @staticmethod
    def eqn_3_3__H_2(H_1: float, P: float, P_P: float, **kwargs):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_2 = H_1 - P + P_P
        result.append(H_2)
        return result

    @staticmethod
    def eqn_3_3__P(H_1: float, H_2: float, P_P: float, **kwargs):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P = H_1 - H_2 + P_P
        result.append(P)
        return result

    @staticmethod
    def eqn_3_3__P_P(H_1: float, H_2: float, P: float, **kwargs):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P_P = -H_1 + H_2 + P
        result.append(P_P)
        return result

