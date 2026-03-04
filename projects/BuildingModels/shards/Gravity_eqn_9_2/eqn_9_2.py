from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_9_2__G_cap import eqn_9_2__G
from .eqn_9_2__M_cap import eqn_9_2__M
from .eqn_9_2__R_cap import eqn_9_2__R
from .eqn_9_2__T_cap import eqn_9_2__T
from .eqn_9_2__pi import eqn_9_2__pi


class Gravity:
    eqn_9_2__G = eqn_9_2__G
    eqn_9_2__M = eqn_9_2__M
    eqn_9_2__R = eqn_9_2__R
    eqn_9_2__T = eqn_9_2__T
    eqn_9_2__pi = eqn_9_2__pi

    @kwasak
    def eqn_9_2(self, G=None, M=None, R=None, T=None, pi=None):
        """
        G := gravitational constant
        M := mass of the Sun
        R := radius of the orbit
        T := orbital period
        """
        return
