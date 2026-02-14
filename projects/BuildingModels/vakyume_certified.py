from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Kinematics:
    @kwasak_static
    def eqn_1_1(a=None, t=None, v=None, v0=None, **kwargs):
        return

    @staticmethod
    def eqn_1_1__a(t: float, v: float, v0: float):
        # [.pyeqn] v = v0 + a * t
        result = []
        a = (v - v0)/t
        result.append(a)
        return result

    @staticmethod
    def eqn_1_1__t(a: float, v: float, v0: float):
        # [.pyeqn] v = v0 + a * t
        result = []
        t = (v - v0)/a
        result.append(t)
        return result

    @staticmethod
    def eqn_1_1__v(a: float, t: float, v0: float):
        # [.pyeqn] v = v0 + a * t
        result = []
        v = a*t + v0
        result.append(v)
        return result

    @staticmethod
    def eqn_1_1__v0(a: float, t: float, v: float):
        # [.pyeqn] v = v0 + a * t
        result = []
        v0 = -a*t + v
        result.append(v0)
        return result


