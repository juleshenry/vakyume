from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Kinematics:
    @kwasak_static
    def eqn_1_3(a=None, dx=None, v=None, v0=None, **kwargs):
        return

    @staticmethod
    def eqn_1_3__a(dx: float, v: float, v0: float):
        # [.pyeqn] v**2 = v0**2 + 2 * a * dx
        result = []
        a = (v**2 - v0**2)/(2*dx)
        result.append(a)
        return result

    @staticmethod
    def eqn_1_3__dx(a: float, v: float, v0: float):
        # [.pyeqn] v**2 = v0**2 + 2 * a * dx
        result = []
        dx = (v**2 - v0**2)/(2*a)
        result.append(dx)
        return result

    @staticmethod
    def eqn_1_3__v(a: float, dx: float, v0: float):
        # [.pyeqn] v**2 = v0**2 + 2 * a * dx
        result = []
        v = -sqrt(2*a*dx + v0**2)
        result.append(v)
        v = sqrt(2*a*dx + v0**2)
        result.append(v)
        return result

    @staticmethod
    def eqn_1_3__v0(a: float, dx: float, v: float):
        # [.pyeqn] v**2 = v0**2 + 2 * a * dx
        result = []
        v0 = -sqrt(-2*a*dx + v**2)
        result.append(v0)
        v0 = sqrt(-2*a*dx + v**2)
        result.append(v0)
        return result

