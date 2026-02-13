from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_06(
        P_1: float = None,
        P_2: float = None,
        S_a: float = None,
        V: float = None,
        t: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_10_06__P_1(P_2: float, S_a: float, V: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        P_1 = P_2 * exp(S_a * t / V)
        result.append(P_1)
        return result

    @staticmethod
    def eqn_10_06__P_2(P_1: float, S_a: float, V: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        P_2 = P_1 * exp(-S_a * t / V)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_10_06__S_a(P_1: float, P_2: float, V: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        S_a = V * log(P_1 / P_2) / t
        result.append(S_a)
        return result

    @staticmethod
    def eqn_10_06__V(P_1: float, P_2: float, S_a: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        V = S_a * t / log(P_1 / P_2)
        result.append(V)
        return result

    @staticmethod
    def eqn_10_06__t(P_1: float, P_2: float, S_a: float, V: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        t = V * log(P_1 / P_2) / S_a
        result.append(t)
        return result


