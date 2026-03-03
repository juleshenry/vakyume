from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class Kinematics:
    @kwasak
    def eqn_1_1(a=None, t=None, v=None, v0=None, **kwargs):
        return
    @staticmethod
    def eqn_1_1__a(t: float, v: float, v0: float):
        # v = v0 + a * t
        result = []
        a = (v - v0) / t
        result.append(a)
        return result
    @staticmethod
    def eqn_1_1__t(a: float, v: float, v0: float):
        # v = v0 + a * t
        result = []
        t = (v - v0) / a
        result.append(t)
        return result
    @staticmethod
    def eqn_1_1__v(a: float, t: float, v0: float):
        # v = v0 + a * t
        result = []
        v = a * t + v0
        result.append(v)
        return result
    @staticmethod
    def eqn_1_1__v0(a: float, t: float, v: float):
        # v = v0 + a * t
        result = []
        v0 = -a * t + v
        result.append(v0)
        return result
    @kwasak
    def eqn_1_2(a=None, t=None, v0=None, x=None, x0=None, **kwargs):
        return
    @staticmethod
    def eqn_1_2__a(t: float, v0: float, x: float, x0: float):
        # x = x0 + v0 * t + 0.5 * a * t**2
        result = []
        a = 2.0 * (-t * v0 + x - x0) / t**2
        result.append(a)
        return result
    @staticmethod
    def eqn_1_2__t(a: float, v0: float, x: float, x0: float):
        # x = x0 + v0 * t + 0.5 * a * t**2
        result = []
        t = (-v0 - sqrt(2.0 * a * x - 2.0 * a * x0 + v0**2)) / a
        result.append(t)
        t = (-v0 + sqrt(2.0 * a * x - 2.0 * a * x0 + v0**2)) / a
        result.append(t)
        return result
    @staticmethod
    def eqn_1_2__v0(a: float, t: float, x: float, x0: float):
        # x = x0 + v0 * t + 0.5 * a * t**2
        result = []
        v0 = (-0.5 * a * t**2 + x - x0) / t
        result.append(v0)
        return result
    @staticmethod
    def eqn_1_2__x(a: float, t: float, v0: float, x0: float):
        # x = x0 + v0 * t + 0.5 * a * t**2
        result = []
        x = 0.5 * a * t**2 + t * v0 + x0
        result.append(x)
        return result
    @staticmethod
    def eqn_1_2__x0(a: float, t: float, v0: float, x: float):
        # x = x0 + v0 * t + 0.5 * a * t**2
        result = []
        x0 = -0.5 * a * t**2 - t * v0 + x
        result.append(x0)
        return result
    @kwasak
    def eqn_1_3(a=None, dx=None, v=None, v0=None, **kwargs):
        return
    @staticmethod
    def eqn_1_3__a(dx: float, v: float, v0: float):
        # v**2 = v0**2 + 2 * a * dx
        result = []
        a = (v**2 - v0**2) / (2 * dx)
        result.append(a)
        return result
    @staticmethod
    def eqn_1_3__dx(a: float, v: float, v0: float):
        # v**2 = v0**2 + 2 * a * dx
        result = []
        dx = (v**2 - v0**2) / (2 * a)
        result.append(dx)
        return result
    @staticmethod
    def eqn_1_3__v(a: float, dx: float, v0: float):
        # v**2 = v0**2 + 2 * a * dx
        result = []
        v = -sqrt(2 * a * dx + v0**2)
        result.append(v)
        v = sqrt(2 * a * dx + v0**2)
        result.append(v)
        return result
    @staticmethod
    def eqn_1_3__v0(a: float, dx: float, v: float):
        # v**2 = v0**2 + 2 * a * dx
        result = []
        v0 = -sqrt(-2 * a * dx + v**2)
        result.append(v0)
        v0 = sqrt(-2 * a * dx + v**2)
        result.append(v0)
        return result
