from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_14(T_c: float = None, T_s: float = None, **kwargs):
        return

    @staticmethod
    def eqn_10_14__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 12
        result = []
        T_c = T_s + 12
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_14__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 12
        result = []
        T_s = T_c - 12
        result.append(T_s)
        return result


