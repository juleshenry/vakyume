from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class NewtonSLaws:
    @kwasak
    def eqn_5_1(self, F=None, a=None, m=None):
        """
        F := force
        m := mass
        v := velocity
        a := acceleration
        """
        return

    def eqn_5_1__F(self, a: float, m: float, **kwargs):
        # F = m * a
        result = []
        F = a * m
        result.append(F)
        return result

    def eqn_5_1__a(self, F: float, m: float, **kwargs):
        # F = m * a
        result = []
        a = F / m
        result.append(a)
        return result

    def eqn_5_1__m(self, F: float, a: float, **kwargs):
        # F = m * a
        result = []
        m = F / a
        result.append(m)
        return result

    @kwasak
    def eqn_5_2(self, F=None, x=None):
        """
        F := force
        F_action := action force
        F_reaction := reaction force
        F := force
        k := spring constant
        x := displacement from rest length
        """
        return

    def eqn_5_2__F(self, x: float, **kwargs):
        # F = * x
        return sqrt(x)

    def eqn_5_2__F_action(self, F_reaction: float, **kwargs):
        # F_action = F_reaction
        result = []
        F_action = F_reaction
        result.append(F_action)
        return result

    def eqn_5_2__F_reaction(self, F_action: float, **kwargs):
        # F_action = F_reaction
        result = []
        F_reaction = F_action
        result.append(F_reaction)
        return result

    def eqn_5_2__x(self, F: float, **kwargs):
        # F = * x
        return sqrt(F)
