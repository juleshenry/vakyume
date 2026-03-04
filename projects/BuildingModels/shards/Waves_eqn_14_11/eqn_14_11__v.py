from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_14_11__v(self, fn: float, n: float, **kwargs):
    # [.pyeqn] fn = n * v / (2L)
    return sqrt(2 * fn * n)
