from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class GaussSLaw:
    @kwasak
    def eqn_17_2(self, E=None, R=None):
        """
        E := electric field
        R := radius
        """
        return

    def eqn_17_2__E(self, R: float, **kwargs):
        # E = E * (4 * R ** 2)
        result = []
        E = 0
        result.append(E)
        return result

    def eqn_17_2__R(self, E: float, **kwargs):
        # E = E * (4 * R ** 2)
        result = []
        R = -1 / 2
        result.append(R)
        R = 1 / 2
        result.append(R)
        return result

    @kwasak
    def eqn_17_3(self, E=None, Q=None, R=None):
        """
        Q := charge
        R := radius
        """
        return

    def eqn_17_3__E(self, Q: float, R: float, **kwargs):
        # E = Q / (4 * R ** 2)
        result = []
        E = Q / (4 * R**2)
        result.append(E)
        return result

    def eqn_17_3__Q(self, E: float, R: float, **kwargs):
        # E = Q / (4 * R ** 2)
        result = []
        Q = 4 * E * R**2
        result.append(Q)
        return result

    def eqn_17_3__R(self, E: float, Q: float, **kwargs):
        # E = Q / (4 * R ** 2)
        result = []
        R = -sqrt(Q / E) / 2
        result.append(R)
        R = sqrt(Q / E) / 2
        result.append(R)
        return result
