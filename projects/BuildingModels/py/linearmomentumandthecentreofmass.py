from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class LinearMomentumAndTheCentreOfMass:
    @kwasak
    def eqn_10_1(self, m=None, p=None, v=None):
        """
        m := mass
        v := velocity
        p := momentum
        t := time
        a := acceleration
        F := force
        """
        return

    def eqn_10_1__m(self, p: float, v: float, **kwargs):
        # p = m * v
        result = []
        m = p / v
        result.append(m)
        return result

    def eqn_10_1__p(self, m: float, v: float, **kwargs):
        # p = m * v
        result = []
        p = m * v
        result.append(p)
        return result

    def eqn_10_1__v(self, m: float, p: float, **kwargs):
        # p = m * v
        result = []
        v = p / m
        result.append(v)
        return result

    @kwasak
    def eqn_10_11(self, F_ext=None, M=None, a_CM=None):
        """
        F_ext := external force
        M := total mass
        a_CM := acceleration of the centre of mass
        """
        return

    def eqn_10_11__F_ext(self, M: float, a_CM: float, **kwargs):
        # F_ext = M * a_CM
        result = []
        F_ext = M * a_CM
        result.append(F_ext)
        return result

    def eqn_10_11__M(self, F_ext: float, a_CM: float, **kwargs):
        # F_ext = M * a_CM
        result = []
        M = F_ext / a_CM
        result.append(M)
        return result

    def eqn_10_11__a_CM(self, F_ext: float, M: float, **kwargs):
        # F_ext = M * a_CM
        result = []
        a_CM = F_ext / M
        result.append(a_CM)
        return result
