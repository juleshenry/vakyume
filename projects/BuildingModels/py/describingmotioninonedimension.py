from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class DescribingMotionInOneDimension:
    @kwasak
    def eqn_3_1(self, t=None, v_x=None, x=None, x_0=None):
        """
        x := position
        x_0 := initial position
        v_x := velocity
        t := time
        """
        return

    def eqn_3_1__t(self, v_x: float, x: float, x_0: float, **kwargs):
        # x = x_0 + v_x * t
        result = []
        t = (x - x_0) / v_x
        result.append(t)
        return result

    def eqn_3_1__v_x(self, t: float, x: float, x_0: float, **kwargs):
        # x = x_0 + v_x * t
        result = []
        v_x = (x - x_0) / t
        result.append(v_x)
        return result

    def eqn_3_1__x(self, t: float, v_x: float, x_0: float, **kwargs):
        # x = x_0 + v_x * t
        result = []
        x = t * v_x + x_0
        result.append(x)
        return result

    def eqn_3_1__x_0(self, t: float, v_x: float, x: float, **kwargs):
        # x = x_0 + v_x * t
        result = []
        x_0 = -t * v_x + x
        result.append(x_0)
        return result

    @kwasak
    def eqn_3_14(self, a=None, a_A=None):
        """
        a := relative acceleration
        a_A := acceleration of the passenger
        """
        return

    def eqn_3_14__a(self, a_A: float, **kwargs):
        # a = a_A
        result = []
        a = a_A
        result.append(a)
        return result

    def eqn_3_14__a_A(self, a: float, **kwargs):
        # a = a_A
        result = []
        a_A = a
        result.append(a_A)
        return result

    @kwasak
    def eqn_3_15(self, t=None, v=None, v_A=None, v_B=None):
        """
        v := velocity of the passenger
        v_B := velocity of the boat
        v_A := velocity of the passenger
        t := time
        """
        return

    def eqn_3_15__t(self, v: float, v_A: float, v_B: float, **kwargs):
        # v = (v_B + v_A) * t
        result = []
        t = v / (v_A + v_B)
        result.append(t)
        return result

    def eqn_3_15__v(self, t: float, v_A: float, v_B: float, **kwargs):
        # v = (v_B + v_A) * t
        result = []
        v = t * (v_A + v_B)
        result.append(v)
        return result

    def eqn_3_15__v_A(self, t: float, v: float, v_B: float, **kwargs):
        # v = (v_B + v_A) * t
        result = []
        v_A = -v_B + v / t
        result.append(v_A)
        return result

    def eqn_3_15__v_B(self, t: float, v: float, v_A: float, **kwargs):
        # v = (v_B + v_A) * t
        result = []
        v_B = -v_A + v / t
        result.append(v_B)
        return result

    @kwasak
    def eqn_3_2(self, ax=None, t=None, v=None, v_0x=None):
        """
        v := velocity
        v_0x := initial velocity
        ax := acceleration
        t := time
        """
        return

    def eqn_3_2__ax(self, t: float, v: float, v_0x: float, **kwargs):
        # v = v_0x + ax * t
        result = []
        ax = (v - v_0x) / t
        result.append(ax)
        return result

    def eqn_3_2__t(self, ax: float, v: float, v_0x: float, **kwargs):
        # v = v_0x + ax * t
        result = []
        t = (v - v_0x) / ax
        result.append(t)
        return result

    def eqn_3_2__v(self, ax: float, t: float, v_0x: float, **kwargs):
        # v = v_0x + ax * t
        result = []
        v = ax * t + v_0x
        result.append(v)
        return result

    def eqn_3_2__v_0x(self, ax: float, t: float, v: float, **kwargs):
        # v = v_0x + ax * t
        result = []
        v_0x = -ax * t + v
        result.append(v_0x)
        return result

    @kwasak
    def eqn_3_3(self, ax=None, t=None, v_0x=None, x=None, x_0=None):
        """
        x := position
        x_0 := initial position
        v_0x := initial velocity
        ax := acceleration
        t := time
        """
        return

    def eqn_3_3__ax(self, t: float, v_0x: float, x: float, x_0: float, **kwargs):
        # x = x_0 + v_0x * t + 0.5 * ax * t ** 2
        result = []
        ax = 2.0 * (-t * v_0x + x - x_0) / t**2
        result.append(ax)
        return result

    def eqn_3_3__t(self, ax: float, v_0x: float, x: float, x_0: float, **kwargs):
        # x = x_0 + v_0x * t + 0.5 * ax * t ** 2
        result = []
        t = (-v_0x - sqrt(2.0 * ax * x - 2.0 * ax * x_0 + v_0x**2)) / ax
        result.append(t)
        t = (-v_0x + sqrt(2.0 * ax * x - 2.0 * ax * x_0 + v_0x**2)) / ax
        result.append(t)
        return result

    def eqn_3_3__v_0x(self, ax: float, t: float, x: float, x_0: float, **kwargs):
        # x = x_0 + v_0x * t + 0.5 * ax * t ** 2
        result = []
        v_0x = (-0.5 * ax * t**2 + x - x_0) / t
        result.append(v_0x)
        return result

    def eqn_3_3__x(self, ax: float, t: float, v_0x: float, x_0: float, **kwargs):
        # x = x_0 + v_0x * t + 0.5 * ax * t ** 2
        result = []
        x = 0.5 * ax * t**2 + t * v_0x + x_0
        result.append(x)
        return result

    def eqn_3_3__x_0(self, ax: float, t: float, v_0x: float, x: float, **kwargs):
        # x = x_0 + v_0x * t + 0.5 * ax * t ** 2
        result = []
        x_0 = -0.5 * ax * t**2 - t * v_0x + x
        result.append(x_0)
        return result

    @kwasak
    def eqn_3_4(self, a=None, t=None, v=None, v_0=None):
        """
        x := position
        t := time
        v_0 := initial velocity
        a := acceleration
        """
        return

    def eqn_3_4__a(self, t: float, v: float, v_0: float, **kwargs):
        # v = v_0 + a * t
        result = []
        a = (v - v_0) / t
        result.append(a)
        return result

    def eqn_3_4__t(self, a: float, v: float, v_0: float, **kwargs):
        # v = v_0 + a * t
        result = []
        t = (v - v_0) / a
        result.append(t)
        return result

    def eqn_3_4__v(self, a: float, t: float, v_0: float, **kwargs):
        # v = v_0 + a * t
        result = []
        v = a * t + v_0
        result.append(v)
        return result

    def eqn_3_4__v_0(self, a: float, t: float, v: float, **kwargs):
        # v = v_0 + a * t
        result = []
        v_0 = -a * t + v
        result.append(v_0)
        return result

    @kwasak
    def eqn_3_6(self, t=None, v_0=None, x=None, x_0=None):
        """
        x := position
        x_0 := initial position
        v_0 := initial velocity
        a := acceleration
        t := time
        """
        return

    def eqn_3_6__t(self, v_0: float, x: float, x_0: float, **kwargs):
        # x = v_0 * t + x_0
        result = []
        t = (x - x_0) / v_0
        result.append(t)
        return result

    def eqn_3_6__v_0(self, t: float, x: float, x_0: float, **kwargs):
        # x = v_0 * t + x_0
        result = []
        v_0 = (x - x_0) / t
        result.append(v_0)
        return result

    def eqn_3_6__x(self, t: float, v_0: float, x_0: float, **kwargs):
        # x = v_0 * t + x_0
        result = []
        x = t * v_0 + x_0
        result.append(x)
        return result

    def eqn_3_6__x_0(self, t: float, v_0: float, x: float, **kwargs):
        # x = v_0 * t + x_0
        result = []
        x_0 = -t * v_0 + x
        result.append(x_0)
        return result
