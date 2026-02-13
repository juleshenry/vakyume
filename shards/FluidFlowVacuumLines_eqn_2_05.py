from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_05(
        D: float = None,
        L: float = None,
        delta_P: float = None,
        mu: float = None,
        q: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_2_05__D(L: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        try:
            D = 2.52647511098426 * (abs(L * mu * q / delta_P)) ** (1 / 4)
            result.append(D)
        except:
            pass
        return result

    @staticmethod
    def eqn_2_05__L(D: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        L = 0.0245436926061703 * D**4 * delta_P / (mu * q)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_05__delta_P(D: float, L: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        delta_P = 40.7436654315252 * L * mu * q / D**4
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_05__mu(D: float, L: float, delta_P: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        mu = 0.0245436926061703 * D**4 * delta_P / (L * q)
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_05__q(D: float, L: float, delta_P: float, mu: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        q = 0.0245436926061703 * D**4 * delta_P / (L * mu)
        result.append(q)
        return result


