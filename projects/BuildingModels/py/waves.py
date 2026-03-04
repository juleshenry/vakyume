from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class Waves:
    @kwasak
    def eqn_14_11(self, fn=None, n=None, v=None):
        """
        n := harmonic number
        L := length of the string
        v := speed of the waves on the string
        """
        return

    def eqn_14_11__fn(self, n: float, v: float, **kwargs):
        # fn = n * v / (2L)
        return [n * v / (2 * sqrt(2))]

    def eqn_14_11__n(self, fn: float, v: float, **kwargs):
        # fn = n * v / (2L)
        return sqrt(2 * v / fn)

    def eqn_14_11__v(self, fn: float, n: float, **kwargs):
        # fn = n * v / (2L)
        return sqrt(2 * fn * n)
