from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class PressMgmt:
    @kwasak_static
    def eqn_2_17(
        L: float = None,
        d: float = None,
        delta_P: float = None,
        mu: float = None,
        q: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_2_17__L(d: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        L = 9.52380952380952 * d**4 * delta_P / (mu * q)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_17__d(L: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        d = 0.569242509762222 * (L * mu * q / delta_P) ** (1 / 4)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_17__delta_P(L: float, d: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        delta_P = 0.105 * L * mu * q / d**4
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_17__mu(L: float, d: float, delta_P: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        mu = 9.52380952380952 * d**4 * delta_P / (L * q)
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_17__q(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        q = 9.52380952380952 * d**4 * delta_P / (L * mu)
        result.append(q)
        return result




