from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class PressMgmt:
    @kwasak_static
    def eqn_3_8(A_C=None, H_2=None, V_P=None, **kwargs):
        return

    @staticmethod
    def eqn_3_8__A_C(H_2: float, V_P: float):
        # [.pyeqn] V_P = A_C * H_2
        result = []
        A_C = V_P/H_2
        result.append(A_C)
        return result

    @staticmethod
    def eqn_3_8__H_2(A_C: float, V_P: float):
        # [.pyeqn] V_P = A_C * H_2
        result = []
        H_2 = V_P/A_C
        result.append(H_2)
        return result

    @staticmethod
    def eqn_3_8__V_P(A_C: float, H_2: float):
        # [.pyeqn] V_P = A_C * H_2
        result = []
        V_P = A_C*H_2
        result.append(V_P)
        return result

