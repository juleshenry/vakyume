from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class Gravity:
    @kwasak
    def eqn_9_2(self, G=None, M=None, R=None, T=None, pi=None):
        """
        G := gravitational constant
        M := mass of the Sun
        R := radius of the orbit
        T := orbital period
        """
        return

    def eqn_9_2__G(self, M: float, R: float, T: float, pi: float, **kwargs):
        # T ** 2 = G * M * R ** 3 / (4 * pi ** 2)
        result = []
        G = 4 * T**2 * pi**2 / (M * R**3)
        result.append(G)
        return result

    def eqn_9_2__M(self, G: float, R: float, T: float, pi: float, **kwargs):
        # T ** 2 = G * M * R ** 3 / (4 * pi ** 2)
        result = []
        M = 4 * T**2 * pi**2 / (G * R**3)
        result.append(M)
        return result

    def eqn_9_2__R(self, G: float, M: float, T: float, pi: float, **kwargs):
        # T ** 2 = G * M * R ** 3 / (4 * pi ** 2)
        result = []
        R = 2 ** (2 / 3) * (T**2 * pi**2 / (G * M)) ** (1 / 3)
        result.append(R)
        R = (
            -(2 ** (2 / 3)) * (T**2 * pi**2 / (G * M)) ** (1 / 3) / 2
            - 2 ** (2 / 3) * sqrt(3) * I * (T**2 * pi**2 / (G * M)) ** (1 / 3) / 2
        )
        result.append(R)
        R = (
            -(2 ** (2 / 3)) * (T**2 * pi**2 / (G * M)) ** (1 / 3) / 2
            + 2 ** (2 / 3) * sqrt(3) * I * (T**2 * pi**2 / (G * M)) ** (1 / 3) / 2
        )
        result.append(R)
        return result

    def eqn_9_2__T(self, G: float, M: float, R: float, pi: float, **kwargs):
        # T ** 2 = G * M * R ** 3 / (4 * pi ** 2)
        result = []
        T = -sqrt(G * M * R**3) / (2 * pi)
        result.append(T)
        T = sqrt(G * M * R**3) / (2 * pi)
        result.append(T)
        return result

    def eqn_9_2__pi(self, G: float, M: float, R: float, T: float, **kwargs):
        # T ** 2 = G * M * R ** 3 / (4 * pi ** 2)
        result = []
        pi = -sqrt(G * M * R**3) / (2 * T)
        result.append(pi)
        pi = sqrt(G * M * R**3) / (2 * T)
        result.append(pi)
        return result

    @kwasak
    def eqn_9_3(self, G=None, M=None, R=None, a=None, m=None):
        """
        m := mass of the object
        a := acceleration due to gravity
        G := gravitational constant
        M := mass of the Earth
        R := radius of the Earth
        """
        return

    def eqn_9_3__G(self, M: float, R: float, a: float, m: float, **kwargs):
        # a = G * M * m / R ** 2
        result = []
        G = R**2 * a / (M * m)
        result.append(G)
        return result

    def eqn_9_3__M(self, G: float, R: float, a: float, m: float, **kwargs):
        # a = G * M * m / R ** 2
        result = []
        M = R**2 * a / (G * m)
        result.append(M)
        return result

    def eqn_9_3__R(self, G: float, M: float, a: float, m: float, **kwargs):
        # a = G * M * m / R ** 2
        result = []
        R = -sqrt(G * M * m / a)
        result.append(R)
        R = sqrt(G * M * m / a)
        result.append(R)
        return result

    def eqn_9_3__a(self, G: float, M: float, R: float, m: float, **kwargs):
        # a = G * M * m / R ** 2
        result = []
        a = G * M * m / R**2
        result.append(a)
        return result

    def eqn_9_3__m(self, G: float, M: float, R: float, a: float, **kwargs):
        # a = G * M * m / R ** 2
        result = []
        m = R**2 * a / (G * M)
        result.append(m)
        return result

    @kwasak
    def eqn_9_4(self, G=None, M=None, U=None, m=None, r=None):
        """
        r := distance
        G := gravitational constant
        M := mass of large body (e.g. Earth)
        m := mass of smaller body (e.g. satellite)
        """
        return

    def eqn_9_4__G(self, M: float, U: float, m: float, r: float, **kwargs):
        # U = - G * M * m / r
        result = []
        G = -U * r / (M * m)
        result.append(G)
        return result

    def eqn_9_4__M(self, G: float, U: float, m: float, r: float, **kwargs):
        # U = - G * M * m / r
        result = []
        M = -U * r / (G * m)
        result.append(M)
        return result

    def eqn_9_4__U(self, G: float, M: float, m: float, r: float, **kwargs):
        # U = - G * M * m / r
        result = []
        U = -G * M * m / r
        result.append(U)
        return result

    def eqn_9_4__m(self, G: float, M: float, U: float, r: float, **kwargs):
        # U = - G * M * m / r
        result = []
        m = -U * r / (G * M)
        result.append(m)
        return result

    def eqn_9_4__r(self, G: float, M: float, U: float, m: float, **kwargs):
        # U = - G * M * m / r
        result = []
        r = -G * M * m / U
        result.append(r)
        return result
