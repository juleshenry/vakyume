from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_12(
        L: float = None,
        d: float = None,
        delta_P: float = None,
        f: float = None,
        g: float = None,
        rho: float = None,
        v: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_2_12__L(d: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        L = 0.464037122969838 * d * delta_P * g / (f * rho * v**2)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_12__d(L: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        d = 2.155 * L * f * rho * v**2 / (delta_P * g)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_12__delta_P(L: float, d: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        delta_P = 2.155 * L * f * rho * v**2 / (d * g)
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_12__f(L: float, d: float, delta_P: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        f = 0.464037122969838 * d * delta_P * g / (L * rho * v**2)
        result.append(f)
        return result

    @staticmethod
    def eqn_2_12__g(L: float, d: float, delta_P: float, f: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        g = 2.155 * L * f * rho * v**2 / (d * delta_P)
        result.append(g)
        return result

    @staticmethod
    def eqn_2_12__rho(L: float, d: float, delta_P: float, f: float, g: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        rho = 0.464037122969838 * d * delta_P * g / (L * f * v**2)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_12__v(L: float, d: float, delta_P: float, f: float, g: float, rho: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        v = 0.681202703290172 * sqrt(d * delta_P * g / (L * f * rho))
        result.append(v)
        return result


