from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_05(
        P_0_1: float = None, P_0_2: float = None, alpha_12: float = None, **kwargs
    ):
        return

    @staticmethod
    def eqn_5_05__P_0_1(P_0_2: float, alpha_12: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_1 = P_0_2 * alpha_12
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_05__P_0_2(P_0_1: float, alpha_12: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_2 = P_0_1 / alpha_12
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_05__alpha_12(P_0_1: float, P_0_2: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        alpha_12 = P_0_1 / P_0_2
        result.append(alpha_12)
        return result


