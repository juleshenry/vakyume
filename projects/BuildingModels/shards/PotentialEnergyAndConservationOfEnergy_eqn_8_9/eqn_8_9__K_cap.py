from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_9__K(self, L: float, **kwargs):
    # [.pyeqn] L = K -
    return sqrt(L + 1) - 1
