from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class PressMgmt:
    @kwasak_static
    def eqn_3_9(A_C=None, H_1=None, H_2=None, P=None, V=None, **kwargs):
        return

    @staticmethod
    def eqn_3_9__A_C(H_1: float, H_2: float, P: float, V: float):
        # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        A_C = P*V/(H_2*(-H_1 + H_2 + P))
        result.append(A_C)
        return result

    @staticmethod
    def eqn_3_9__H_1(A_C: float, H_2: float, P: float, V: float):
        # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        H_1 = H_2 + P - P*V/(A_C*H_2)
        result.append(H_1)
        return result

    @staticmethod
    def eqn_3_9__H_2(A_C: float, H_1: float, P: float, V: float):
        # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        H_2 = (A_C*(H_1 - P) - sqrt(A_C*(A_C*H_1**2 - 2*A_C*H_1*P + A_C*P**2 + 4*P*V)))/(2*A_C)
        result.append(H_2)
        H_2 = (A_C*(H_1 - P) + sqrt(A_C*(A_C*H_1**2 - 2*A_C*H_1*P + A_C*P**2 + 4*P*V)))/(2*A_C)
        result.append(H_2)
        return result

    @staticmethod
    def eqn_3_9__P(A_C: float, H_1: float, H_2: float, V: float):
        # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        P = A_C*H_2*(H_1 - H_2)/(A_C*H_2 - V)
        result.append(P)
        return result

    @staticmethod
    def eqn_3_9__V(A_C: float, H_1: float, H_2: float, P: float):
        # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        V = A_C*H_2*(-H_1 + H_2 + P)/P
        result.append(V)
        return result

