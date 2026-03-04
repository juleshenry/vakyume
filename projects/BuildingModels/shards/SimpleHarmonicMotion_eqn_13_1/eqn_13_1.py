from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_13_1__E_cap import eqn_13_1__E
from .eqn_13_1__m import eqn_13_1__m
from .eqn_13_1__v import eqn_13_1__v
from .eqn_13_1__x import eqn_13_1__x


class SimpleHarmonicMotion:
    eqn_13_1__E = eqn_13_1__E
    eqn_13_1__m = eqn_13_1__m
    eqn_13_1__v = eqn_13_1__v
    eqn_13_1__x = eqn_13_1__x

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
