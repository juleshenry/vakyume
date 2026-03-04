from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_9_4__G_cap import eqn_9_4__G
from .eqn_9_4__M_cap import eqn_9_4__M
from .eqn_9_4__U_cap import eqn_9_4__U
from .eqn_9_4__m import eqn_9_4__m
from .eqn_9_4__r import eqn_9_4__r


class Gravity:
    eqn_9_4__G = eqn_9_4__G
    eqn_9_4__M = eqn_9_4__M
    eqn_9_4__U = eqn_9_4__U
    eqn_9_4__m = eqn_9_4__m
    eqn_9_4__r = eqn_9_4__r

    @kwasak
    def eqn_9_4(self, G=None, M=None, U=None, m=None, r=None):
        """
        r := distance
        G := gravitational constant
        M := mass of large body (e.g. Earth)
        m := mass of smaller body (e.g. satellite)
        """
        return
