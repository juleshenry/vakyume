from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class DescribingMotionInMultipleDimensions:
    @kwasak
    def eqn_4_1(self, t=None, v=None, x_1=None, x_2=None, y_1=None, y_2=None):
        """
        t := time
        x_1 := initial x-coordinate
        y_1 := initial y-coordinate
        x_2 := final x-coordinate
        y_2 := final y-coordinate
        """
        return

    def eqn_4_1__t(
        self, v: float, x_1: float, x_2: float, y_1: float, y_2: float, **kwargs
    ):
        # v = (x_2 - x_1) / t + (y_2 - y_1) / t
        result = []
        t = (-x_1 + x_2 - y_1 + y_2) / v
        result.append(t)
        return result

    def eqn_4_1__v(
        self, t: float, x_1: float, x_2: float, y_1: float, y_2: float, **kwargs
    ):
        # v = (x_2 - x_1) / t + (y_2 - y_1) / t
        result = []
        v = (-x_1 + x_2 - y_1 + y_2) / t
        result.append(v)
        return result

    def eqn_4_1__x_1(
        self, t: float, v: float, x_2: float, y_1: float, y_2: float, **kwargs
    ):
        # v = (x_2 - x_1) / t + (y_2 - y_1) / t
        result = []
        x_1 = -t * v + x_2 - y_1 + y_2
        result.append(x_1)
        return result

    def eqn_4_1__x_2(
        self, t: float, v: float, x_1: float, y_1: float, y_2: float, **kwargs
    ):
        # v = (x_2 - x_1) / t + (y_2 - y_1) / t
        result = []
        x_2 = t * v + x_1 + y_1 - y_2
        result.append(x_2)
        return result

    def eqn_4_1__y_1(
        self, t: float, v: float, x_1: float, x_2: float, y_2: float, **kwargs
    ):
        # v = (x_2 - x_1) / t + (y_2 - y_1) / t
        result = []
        y_1 = -t * v - x_1 + x_2 + y_2
        result.append(y_1)
        return result

    def eqn_4_1__y_2(
        self, t: float, v: float, x_1: float, x_2: float, y_1: float, **kwargs
    ):
        # v = (x_2 - x_1) / t + (y_2 - y_1) / t
        result = []
        y_2 = t * v + x_1 - x_2 + y_1
        result.append(y_2)
        return result

    @kwasak
    def eqn_4_3(self, v=None, vx=None, vy=None, x=None, y=None):
        """
        v := velocity
        x := x-component
        y := y-component
        t := time
        """
        return

    def eqn_4_3__v(self, vx: float, vy: float, x: float, y: float, **kwargs):
        # v = vx * x + vy * y
        result = []
        v = vx * x + vy * y
        result.append(v)
        return result

    def eqn_4_3__vx(self, v: float, vy: float, x: float, y: float, **kwargs):
        # v = vx * x + vy * y
        result = []
        vx = (v - vy * y) / x
        result.append(vx)
        return result

    def eqn_4_3__vy(self, v: float, vx: float, x: float, y: float, **kwargs):
        # v = vx * x + vy * y
        result = []
        vy = (v - vx * x) / y
        result.append(vy)
        return result

    def eqn_4_3__x(self, v: float, vx: float, vy: float, y: float, **kwargs):
        # v = vx * x + vy * y
        result = []
        x = (v - vy * y) / vx
        result.append(x)
        return result

    def eqn_4_3__y(self, v: float, vx: float, vy: float, x: float, **kwargs):
        # v = vx * x + vy * y
        result = []
        y = (v - vx * x) / vy
        result.append(y)
        return result

    @kwasak
    def eqn_4_4(self, a=None, ax=None, ay=None, x=None, y=None):
        """
        a := acceleration
        x := x-component
        y := y-component
        t := time
        """
        return

    def eqn_4_4__a(self, ax: float, ay: float, x: float, y: float, **kwargs):
        # a = ax * x + ay * y
        result = []
        a = ax * x + ay * y
        result.append(a)
        return result

    def eqn_4_4__ax(self, a: float, ay: float, x: float, y: float, **kwargs):
        # a = ax * x + ay * y
        result = []
        ax = (a - ay * y) / x
        result.append(ax)
        return result

    def eqn_4_4__ay(self, a: float, ax: float, x: float, y: float, **kwargs):
        # a = ax * x + ay * y
        result = []
        ay = (a - ax * x) / y
        result.append(ay)
        return result

    def eqn_4_4__x(self, a: float, ax: float, ay: float, y: float, **kwargs):
        # a = ax * x + ay * y
        result = []
        x = (a - ay * y) / ax
        result.append(x)
        return result

    def eqn_4_4__y(self, a: float, ax: float, ay: float, x: float, **kwargs):
        # a = ax * x + ay * y
        result = []
        y = (a - ax * x) / ay
        result.append(y)
        return result

    @kwasak
    def eqn_4_5(self, a=None, t=None, v=None, v_0=None):
        """
        v := velocity
        v_0 := initial velocity
        a := acceleration
        t := time
        """
        return

    def eqn_4_5__a(self, t: float, v: float, v_0: float, **kwargs):
        # v = v_0 + a * t
        result = []
        a = (v - v_0) / t
        result.append(a)
        return result

    def eqn_4_5__t(self, a: float, v: float, v_0: float, **kwargs):
        # v = v_0 + a * t
        result = []
        t = (v - v_0) / a
        result.append(t)
        return result

    def eqn_4_5__v(self, a: float, t: float, v_0: float, **kwargs):
        # v = v_0 + a * t
        result = []
        v = a * t + v_0
        result.append(v)
        return result

    def eqn_4_5__v_0(self, a: float, t: float, v: float, **kwargs):
        # v = v_0 + a * t
        result = []
        v_0 = -a * t + v
        result.append(v_0)
        return result

    @kwasak
    def eqn_4_6(self, a=None, r=None, r_0=None, t=None, v_0=None):
        """
        r := position
        r_0 := initial position
        v_0 := initial velocity
        a := acceleration
        t := time
        """
        return

    def eqn_4_6__a(self, r: float, r_0: float, t: float, v_0: float, **kwargs):
        # r = r_0 + v_0 * t + 0.5 * a * t ** 2
        result = []
        a = 2.0 * (r - r_0 - t * v_0) / t**2
        result.append(a)
        return result

    def eqn_4_6__r(self, a: float, r_0: float, t: float, v_0: float, **kwargs):
        # r = r_0 + v_0 * t + 0.5 * a * t ** 2
        result = []
        r = 0.5 * a * t**2 + r_0 + t * v_0
        result.append(r)
        return result

    def eqn_4_6__r_0(self, a: float, r: float, t: float, v_0: float, **kwargs):
        # r = r_0 + v_0 * t + 0.5 * a * t ** 2
        result = []
        r_0 = -0.5 * a * t**2 + r - t * v_0
        result.append(r_0)
        return result

    def eqn_4_6__t(self, a: float, r: float, r_0: float, v_0: float, **kwargs):
        # r = r_0 + v_0 * t + 0.5 * a * t ** 2
        result = []
        t = (-v_0 - sqrt(2.0 * a * r - 2.0 * a * r_0 + v_0**2)) / a
        result.append(t)
        t = (-v_0 + sqrt(2.0 * a * r - 2.0 * a * r_0 + v_0**2)) / a
        result.append(t)
        return result

    def eqn_4_6__v_0(self, a: float, r: float, r_0: float, t: float, **kwargs):
        # r = r_0 + v_0 * t + 0.5 * a * t ** 2
        result = []
        v_0 = (-0.5 * a * t**2 + r - r_0) / t
        result.append(v_0)
        return result
