from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_14_11__fn import eqn_14_11__fn
from .eqn_14_11__n import eqn_14_11__n
from .eqn_14_11__v import eqn_14_11__v


class Waves:
    eqn_14_11__fn = eqn_14_11__fn
    eqn_14_11__n = eqn_14_11__n
    eqn_14_11__v = eqn_14_11__v

    @kwasak
    def eqn_14_11(self, fn=None, n=None, v=None):
        """
        n := harmonic number
        L := length of the string
        v := speed of the waves on the string
        """
        return
