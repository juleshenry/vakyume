from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class RotaryPistonVane:
    @kwasak_static
    def eqn_11_3(F_s=None, t=None, t_c=None, **kwargs):
        return

    @staticmethod
    def eqn_11_3__F_s(t: float, t_c: float, **kwargs):
        # [.pyeqn] t = t_c * F_s
        result = []
        F_s = t/t_c
        result.append(F_s)
        return result

    @staticmethod
    def eqn_11_3__t(F_s: float, t_c: float, **kwargs):
        # [.pyeqn] t = t_c * F_s
        result = []
        t = F_s*t_c
        result.append(t)
        return result

    @staticmethod
    def eqn_11_3__t_c(F_s: float, t: float, **kwargs):
        # [.pyeqn] t = t_c * F_s
        result = []
        t_c = t/F_s
        result.append(t_c)
        return result

