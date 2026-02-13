from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class VacuumTheory:
    @kwasak_static
    def eqn_1_08(
        M: float = None,
        P: float = None,
        R: float = None,
        T: float = None,
        V: float = None,
        m: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_1_08__M(P: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        M = R * T * m / (P * V)
        result.append(M)
        return result

    @staticmethod
    def eqn_1_08__P(M: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        P = R * T * m / (M * V)
        result.append(P)
        return result

    @staticmethod
    def eqn_1_08__R(M: float, P: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        R = M * P * V / (T * m)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_08__T(M: float, P: float, R: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        T = M * P * V / (R * m)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_08__V(M: float, P: float, R: float, T: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        V = R * T * m / (M * P)
        result.append(V)
        return result

    @staticmethod
    def eqn_1_08__m(M: float, P: float, R: float, T: float, V: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        m = M * P * V / (R * T)
        result.append(m)
        return result


