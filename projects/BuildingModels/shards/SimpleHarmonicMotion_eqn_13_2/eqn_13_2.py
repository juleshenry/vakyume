from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_13_2__x0 import eqn_13_2__x0
from .eqn_13_2__x1 import eqn_13_2__x1
from .eqn_13_2__x2 import eqn_13_2__x2


class SimpleHarmonicMotion:
    eqn_13_2__x0 = eqn_13_2__x0
    eqn_13_2__x1 = eqn_13_2__x1
    eqn_13_2__x2 = eqn_13_2__x2

    @kwasak
    def eqn_13_2(self, x0=None, x1=None, x2=None):
        """
        m := mass
        k := spring constant
        x := position
        t := time
        m := mass
        k1 := spring constant 1
        k2 := spring constant 2
        x0 := equilibrium position
        x1 := position of spring 1 at rest
        x2 := position of spring 2 at rest
        """
        return
