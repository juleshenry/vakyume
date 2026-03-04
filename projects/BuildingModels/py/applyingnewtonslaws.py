from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class ApplyingNewtonSLaws:
    @kwasak
    def eqn_6_1(self, F=None, F_N=None):
        """
        F := force
        F_N := normal force
        """
        return

    def eqn_6_1__F(self, F_N: float, **kwargs):
        # F = F_N
        result = []
        F = F_N
        result.append(F)
        return result

    def eqn_6_1__F_N(self, F: float, **kwargs):
        # F = F_N
        result = []
        F_N = F
        result.append(F_N)
        return result

    @kwasak
    def eqn_6_2(self, a_2=None, m=None, x=None):
        """
        x := horizontal position
        x_0 := initial horizontal position
        a_2 := acceleration
        m := mass
        k_2 := kinetic friction coefficient
        N_2 := normal force
        g := gravitational acceleration
        """
        return

    def eqn_6_2__a_2(self, m: float, x: float, **kwargs):
        # x = - m * a_2
        result = []
        a_2 = -x / m
        result.append(a_2)
        return result

    def eqn_6_2__m(self, a_2: float, x: float, **kwargs):
        # x = - m * a_2
        result = []
        m = -x / a_2
        result.append(m)
        return result

    def eqn_6_2__x(self, a_2: float, m: float, **kwargs):
        # x = - m * a_2
        result = []
        x = -a_2 * m
        result.append(x)
        return result
