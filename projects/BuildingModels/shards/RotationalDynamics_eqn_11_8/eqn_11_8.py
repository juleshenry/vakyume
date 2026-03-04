from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_11_8__I_cap import eqn_11_8__I
from .eqn_11_8__i import eqn_11_8__i
from .eqn_11_8__m import eqn_11_8__m
from .eqn_11_8__r import eqn_11_8__r


class RotationalDynamics:
    eqn_11_8__I = eqn_11_8__I
    eqn_11_8__i = eqn_11_8__i
    eqn_11_8__m = eqn_11_8__m
    eqn_11_8__r = eqn_11_8__r

    @kwasak
    def eqn_11_8(self, I=None, i=None, m=None, r=None):
        """
        I := moment of inertia
        i := moment of inertia
        m := mass
        r := distance
        """
        return
