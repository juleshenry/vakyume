from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_02a(
        K_1: float = None, K_2: float = None, alpha_1_2: float = None, **kwargs
    ):
        return

    @staticmethod
    def eqn_5_02a__K_1(K_2: float, alpha_1_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        K_1 = K_2 * alpha_1_2
        result.append(K_1)
        return result

    @staticmethod
    def eqn_5_02a__K_2(K_1: float, alpha_1_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        K_2 = K_1 / alpha_1_2
        result.append(K_2)
        return result

    @staticmethod
    def eqn_5_02a__alpha_1_2(K_1: float, K_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        alpha_1_2 = K_1 / K_2
        result.append(alpha_1_2)
        return result


