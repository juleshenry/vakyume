from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_14_11__fn(self, n: float, v: float, **kwargs):
    # [.pyeqn] fn = n * v / (2L)
    return [n * v / (2 * sqrt(2))]
