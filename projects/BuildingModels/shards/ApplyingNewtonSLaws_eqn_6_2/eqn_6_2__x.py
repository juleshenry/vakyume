from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_6_2__x(self, a_2: float, m: float, **kwargs):
    # [.pyeqn] x = - m * a_2
    result = []
    x = -a_2 * m
    result.append(x)
    return result
