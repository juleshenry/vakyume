from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_09(
        T_c: float = None, T_s: float = None, delta_T: float = None, **kwargs
    ):
        return

    @staticmethod
    def eqn_10_09__T_c(T_s: float, delta_T: float):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        T_c = T_s + delta_T
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_09__T_s(T_c: float, delta_T: float):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        T_s = T_c - delta_T
        result.append(T_s)
        return result

    @staticmethod
    def eqn_10_09__delta_T(T_c: float, T_s: float):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        delta_T = T_c - T_s
        result.append(delta_T)
        return result


