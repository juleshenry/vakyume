from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_07(
        P: float = None,
        P_0_i: float = None,
        gamma_i: float = None,
        x_i: float = None,
        y_i: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_5_07__P(P_0_i: float, gamma_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P = P_0_i * gamma_i * x_i / y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_5_07__P_0_i(P: float, gamma_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P_0_i = P * y_i / (gamma_i * x_i)
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_07__gamma_i(P: float, P_0_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        gamma_i = P * y_i / (P_0_i * x_i)
        result.append(gamma_i)
        return result

    @staticmethod
    def eqn_5_07__x_i(P: float, P_0_i: float, gamma_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        x_i = P * y_i / (P_0_i * gamma_i)
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_07__y_i(P: float, P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        y_i = P_0_i * gamma_i * x_i / P
        result.append(y_i)
        return result


