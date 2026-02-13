from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SteamJetInjectors:
    @kwasak_static
    def eqn_8_09(
        E_j: float = None,
        E_m: float = None,
        e: float = None,
        r: float = None,
        s: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_8_09__E_j(E_m: float, e: float, r: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        E_j = 0.341296928327645 * E_m * r * s / e
        result.append(E_j)
        return result

    @staticmethod
    def eqn_8_09__E_m(E_j: float, e: float, r: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        E_m = 2.93 * E_j * e / (r * s)
        result.append(E_m)
        return result

    @staticmethod
    def eqn_8_09__e(E_j: float, E_m: float, r: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        e = 0.341296928327645 * E_m * r * s / E_j
        result.append(e)
        return result

    @staticmethod
    def eqn_8_09__r(E_j: float, E_m: float, e: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        r = 2.93 * E_j * e / (E_m * s)
        result.append(r)
        return result

    @staticmethod
    def eqn_8_09__s(E_j: float, E_m: float, e: float, r: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        s = 2.93 * E_j * e / (E_m * r)
        result.append(s)
        return result




