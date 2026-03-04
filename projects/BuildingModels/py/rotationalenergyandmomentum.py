from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class RotationalEnergyAndMomentum:
    @kwasak
    def eqn_12_11(self, I=None, L=None, v=None):
        """
        m := mass of a particle
        r := position of the particle
        v := velocity of the particle
        I := moment of inertia
        """
        return

    def eqn_12_11__I(self, L: float, v: float, **kwargs):
        # L = (1 / 2) * I * v ** 2
        result = []
        I = 2 * L / v**2
        result.append(I)
        return result

    def eqn_12_11__L(self, I: float, v: float, **kwargs):
        # L = (1 / 2) * I * v ** 2
        result = []
        L = I * v**2 / 2
        result.append(L)
        return result

    def eqn_12_11__v(self, I: float, L: float, **kwargs):
        # L = (1 / 2) * I * v ** 2
        result = []
        v = -sqrt(2) * sqrt(L / I)
        result.append(v)
        v = sqrt(2) * sqrt(L / I)
        result.append(v)
        return result
