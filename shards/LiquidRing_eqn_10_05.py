from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_05(
        P_1: float = None,
        P_2: float = None,
        S_p: float = None,
        V: float = None,
        t: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_10_05__P_1(P_2: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        P_1 = P_2 * exp(S_p * t / V)
        result.append(P_1)
        return result

    @staticmethod
    def eqn_10_05__P_2(P_1: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        P_2 = P_1 * exp(-S_p * t / V)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_10_05__S_p(P_1: float, P_2: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        S_p = V * log(P_1 / P_2) / t
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_05__V(P_1: float, P_2: float, S_p: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        V = S_p * t / log(P_1 / P_2)
        result.append(V)
        return result

    @staticmethod
    def eqn_10_05__t(P_1: float, P_2: float, S_p: float, V: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        t = V * log(P_1 / P_2) / S_p
        result.append(t)
        return result


