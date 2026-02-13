from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class PressMgmt:
    @kwasak_static
    def eqn_3_05(
        P: float = None, P_P: float = None, V: float = None, V_P: float = None, **kwargs
    ):
        return

    @staticmethod
    def eqn_3_05__P(P_P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P = P_P * V_P / V
        result.append(P)
        return result

    @staticmethod
    def eqn_3_05__P_P(P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P_P = P * V / V_P
        result.append(P_P)
        return result

    @staticmethod
    def eqn_3_05__V(P: float, P_P: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V = P_P * V_P / P
        result.append(V)
        return result

    @staticmethod
    def eqn_3_05__V_P(P: float, P_P: float, V: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V_P = P * V / P_P
        result.append(V_P)
        return result


