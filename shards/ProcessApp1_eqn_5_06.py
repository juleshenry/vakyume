from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_06(
        P_0_i: float = None,
        gamma_i: float = None,
        p_i: float = None,
        x_i: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_5_06__P_0_i(gamma_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        P_0_i = p_i / (gamma_i * x_i)
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_06__gamma_i(P_0_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        gamma_i = p_i / (P_0_i * x_i)
        result.append(gamma_i)
        return result

    @staticmethod
    def eqn_5_06__p_i(P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        p_i = P_0_i * gamma_i * x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_06__x_i(P_0_i: float, gamma_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        x_i = p_i / (P_0_i * gamma_i)
        result.append(x_i)
        return result


