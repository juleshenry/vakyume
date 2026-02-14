from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_12(T_c=None, T_s=None, **kwargs):
        return

    @staticmethod
    def eqn_10_12__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 5
        result = []
        T_c = T_s + 5
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_12__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 5
        result = []
        T_s = T_c - 5
        result.append(T_s)
        return result

