from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_5_1__F_cap import eqn_5_1__F
from .eqn_5_1__a import eqn_5_1__a
from .eqn_5_1__m import eqn_5_1__m


class NewtonSLaws:
    eqn_5_1__F = eqn_5_1__F
    eqn_5_1__a = eqn_5_1__a
    eqn_5_1__m = eqn_5_1__m

    @kwasak
    def eqn_5_1(self, F=None, a=None, m=None):
        """
        F := force
        m := mass
        v := velocity
        a := acceleration
        """
        return
