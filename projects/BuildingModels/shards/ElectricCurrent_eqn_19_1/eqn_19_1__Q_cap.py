from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_19_1__Q(self, I: float, t: float, **kwargs):
    # [.pyeqn] I = Q / t
    result = []
    Q = I * t
    result.append(Q)
    return result
