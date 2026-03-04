from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class RotationalDynamics:
    @kwasak
    def eqn_11_10(self, ICM=None, Ih=None, M=None, h=None):
        """
        Ih := moment of inertia about a parallel axis
        ICM := moment of inertia about the centre of mass
        M := mass of the object
        h := distance from the centre of mass to the parallel axis
        """
        return

    def eqn_11_10__ICM(self, Ih: float, M: float, h: float, **kwargs):
        # Ih = ICM + M * h ** 2
        result = []
        ICM = Ih - M * h**2
        result.append(ICM)
        return result

    def eqn_11_10__Ih(self, ICM: float, M: float, h: float, **kwargs):
        # Ih = ICM + M * h ** 2
        result = []
        Ih = ICM + M * h**2
        result.append(Ih)
        return result

    def eqn_11_10__M(self, ICM: float, Ih: float, h: float, **kwargs):
        # Ih = ICM + M * h ** 2
        result = []
        M = (-ICM + Ih) / h**2
        result.append(M)
        return result

    def eqn_11_10__h(self, ICM: float, Ih: float, M: float, **kwargs):
        # Ih = ICM + M * h ** 2
        result = []
        h = sqrt((-ICM + Ih) / M)
        result.append(h)
        h = -sqrt(-(ICM - Ih) / M)
        result.append(h)
        return result

    @kwasak
    def eqn_11_8(self, I=None, i=None, m=None, r=None):
        """
        I := moment of inertia
        i := moment of inertia
        m := mass
        r := distance
        """
        return

    def eqn_11_8__I(self, i: float, m: float, r: float, **kwargs):
        # I = (m * r ** 2) / i
        result = []
        I = m * r**2 / i
        result.append(I)
        return result

    def eqn_11_8__i(self, I: float, m: float, r: float, **kwargs):
        # I = (m * r ** 2) / i
        result = []
        i = m * r**2 / I
        result.append(i)
        return result

    def eqn_11_8__m(self, I: float, i: float, r: float, **kwargs):
        # I = (m * r ** 2) / i
        result = []
        m = I * i / r**2
        result.append(m)
        return result

    def eqn_11_8__r(self, I: float, i: float, m: float, **kwargs):
        # I = (m * r ** 2) / i
        result = []
        r = -sqrt(I * i / m)
        result.append(r)
        r = sqrt(I * i / m)
        result.append(r)
        return result
