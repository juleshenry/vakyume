from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_6_2__a_2 import eqn_6_2__a_2
from .eqn_6_2__m import eqn_6_2__m
from .eqn_6_2__x import eqn_6_2__x


class ApplyingNewtonSLaws:
    eqn_6_2__a_2 = eqn_6_2__a_2
    eqn_6_2__m = eqn_6_2__m
    eqn_6_2__x = eqn_6_2__x

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
