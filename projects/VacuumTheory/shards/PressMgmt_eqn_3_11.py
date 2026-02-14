from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class PressMgmt:
    @kwasak_static
    def eqn_3_11(A_C=None, H_2=None, P=None, V=None, **kwargs):
        return

    @staticmethod
    def eqn_3_11__A_C(H_2: float, P: float, V: float, **kwargs):
        # [.pyeqn] P = A_C / V * (H_2) ** 2
        result = []
        A_C = P*V/H_2**2
        result.append(A_C)
        return result

    @staticmethod
    def eqn_3_11__H_2(A_C: float, P: float, V: float, **kwargs):
        # [.pyeqn] P = A_C / V * (H_2) ** 2
        result = []
        H_2 = -sqrt(P*V/A_C)
        result.append(H_2)
        H_2 = sqrt(P*V/A_C)
        result.append(H_2)
        return result

    @staticmethod
    def eqn_3_11__P(A_C: float, H_2: float, V: float, **kwargs):
        # [.pyeqn] P = A_C / V * (H_2) ** 2
        result = []
        P = A_C*H_2**2/V
        result.append(P)
        return result

    @staticmethod
    def eqn_3_11__V(A_C: float, H_2: float, P: float, **kwargs):
        # [.pyeqn] P = A_C / V * (H_2) ** 2
        result = []
        V = A_C*H_2**2/P
        result.append(V)
        return result

