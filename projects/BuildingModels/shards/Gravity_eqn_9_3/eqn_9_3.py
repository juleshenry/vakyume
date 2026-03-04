from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_9_3__G_cap import eqn_9_3__G
from .eqn_9_3__M_cap import eqn_9_3__M
from .eqn_9_3__R_cap import eqn_9_3__R
from .eqn_9_3__a import eqn_9_3__a
from .eqn_9_3__m import eqn_9_3__m


class Gravity:
    eqn_9_3__G = eqn_9_3__G
    eqn_9_3__M = eqn_9_3__M
    eqn_9_3__R = eqn_9_3__R
    eqn_9_3__a = eqn_9_3__a
    eqn_9_3__m = eqn_9_3__m

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
