from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class SimpleHarmonicMotion:
    @kwasak
    def eqn_13_1(self, E=None, m=None, v=None, x=None):
        """
        x := position
        A := amplitude
        k := spring constant
        m := mass
        v := velocity
        """
        return

    def eqn_13_1__E(self, m: float, v: float, x: float, **kwargs):
        # E = 0.5 ** x ** 2 + 0.5 * m * v ** 2
        result = []
        E = 0.5 * m * v**2 + 2.0 ** (-(x**2))
        result.append(E)
        return result

    def eqn_13_1__m(self, E: float, v: float, x: float, **kwargs):
        # E = 0.5 ** x ** 2 + 0.5 * m * v ** 2
        result = []
        m = 2.0 * E / v**2 - 2.0 / (2.0 ** (x**2) * v**2)
        result.append(m)
        return result

    def eqn_13_1__v(self, E: float, m: float, x: float, **kwargs):
        # E = 0.5 ** x ** 2 + 0.5 * m * v ** 2
        result = []
        v = -1.4142135623731 * sqrt((E - 1 / 2.0 ** (x**2)) / m)
        result.append(v)
        v = 1.4142135623731 * sqrt((E - 1 / 2.0 ** (x**2)) / m)
        result.append(v)
        return result

    def eqn_13_1__x(self, E: float, m: float, v: float, **kwargs):
        # E = 0.5 ** x ** 2 + 0.5 * m * v ** 2
        result = []
        x = -1.20112240878645 * sqrt(log(2.0 / (2.0 * E - m * v**2)))
        result.append(x)
        x = 1.20112240878645 * sqrt(log(2.0 / (2.0 * E - m * v**2)))
        result.append(x)
        return result

    @kwasak
    def eqn_13_2(self, x0=None, x1=None, x2=None):
        """
        m := mass
        k := spring constant
        x := position
        t := time
        m := mass
        k1 := spring constant 1
        k2 := spring constant 2
        x0 := equilibrium position
        x1 := position of spring 1 at rest
        x2 := position of spring 2 at rest
        """
        return

    def eqn_13_2__d2x(self, dt2: float, m: float, x: float, **kwargs):
        # m * d2x / dt2 = - * x
        # Placeholder for numerical solver
        raise UnsolvedException("Pending LLM/Manual Repair")

    def eqn_13_2__dt2(self, d2x: float, m: float, x: float, **kwargs):
        # m * d2x / dt2 = - * x
        # Placeholder for numerical solver
        raise UnsolvedException("Pending LLM/Manual Repair")

    def eqn_13_2__m(self, d2x: float, dt2: float, x: float, **kwargs):
        # m * d2x / dt2 = - * x
        # Placeholder for numerical solver
        raise UnsolvedException("Pending LLM/Manual Repair")

    def eqn_13_2__x(self, d2x: float, dt2: float, m: float, **kwargs):
        # m * d2x / dt2 = - * x
        # Placeholder for numerical solver
        raise UnsolvedException("Pending LLM/Manual Repair")

    def eqn_13_2__x0(self, x1: float, x2: float, **kwargs):
        # 1 * x1 + 2 * x2 = (1 + 2) * x0
        result = []
        x0 = x1 / 3 + 2 * x2 / 3
        result.append(x0)
        return result

    def eqn_13_2__x1(self, x0: float, x2: float, **kwargs):
        # 1 * x1 + 2 * x2 = (1 + 2) * x0
        result = []
        x1 = 3 * x0 - 2 * x2
        result.append(x1)
        return result

    def eqn_13_2__x2(self, x0: float, x1: float, **kwargs):
        # 1 * x1 + 2 * x2 = (1 + 2) * x0
        result = []
        x2 = 3 * x0 / 2 - x1 / 2
        result.append(x2)
        return result
