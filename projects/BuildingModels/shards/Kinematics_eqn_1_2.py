from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak
import numpy as np


class Kinematics:
    @kwasak
    def eqn_1_2(a=None, t=None, v0=None, x=None, x0=None, **kwargs):
        return

    @staticmethod
    def eqn_1_2__a(t: float, v0: float, x: float, x0: float):
        # [.pyeqn] x = x0 + v0 * t + 0.5 * a * t**2
        result = []
        a = 2.0 * (-t * v0 + x - x0) / t**2
        result.append(a)
        return result

    @staticmethod
    def eqn_1_2__t(a: float, v0: float, x: float, x0: float):
        # [.pyeqn] x = x0 + v0 * t + 0.5 * a * t**2
        result = []
        t = (-v0 - sqrt(2.0 * a * x - 2.0 * a * x0 + v0**2)) / a
        result.append(t)
        t = (-v0 + sqrt(2.0 * a * x - 2.0 * a * x0 + v0**2)) / a
        result.append(t)
        return result

    @staticmethod
    def eqn_1_2__v0(a: float, t: float, x: float, x0: float):
        # [.pyeqn] x = x0 + v0 * t + 0.5 * a * t**2
        result = []
        v0 = (-0.5 * a * t**2 + x - x0) / t
        result.append(v0)
        return result

    @staticmethod
    def eqn_1_2__x(a: float, t: float, v0: float, x0: float):
        # [.pyeqn] x = x0 + v0 * t + 0.5 * a * t**2
        result = []
        x = 0.5 * a * t**2 + t * v0 + x0
        result.append(x)
        return result

    @staticmethod
    def eqn_1_2__x0(a: float, t: float, v0: float, x: float):
        # [.pyeqn] x = x0 + v0 * t + 0.5 * a * t**2
        result = []
        x0 = -0.5 * a * t**2 - t * v0 + x
        result.append(x0)
        return result
