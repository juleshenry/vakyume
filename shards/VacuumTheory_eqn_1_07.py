from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class VacuumTheory:
    @kwasak_static
    def eqn_1_07(
        R: float = None,
        T: float = None,
        V: float = None,
        n: float = None,
        p: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_1_07__R(T: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        R = V * p / (T * n)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_07__T(R: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        T = V * p / (R * n)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_07__V(R: float, T: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        V = R * T * n / p
        result.append(V)
        return result

    @staticmethod
    def eqn_1_07__n(R: float, T: float, V: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        n = V * p / (R * T)
        result.append(n)
        return result

    @staticmethod
    def eqn_1_07__p(R: float, T: float, V: float, n: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        p = R * T * n / V
        result.append(p)
        return result


