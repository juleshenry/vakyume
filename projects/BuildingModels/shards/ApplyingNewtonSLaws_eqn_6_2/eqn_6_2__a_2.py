from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_6_2__a_2(self, m: float, x: float, **kwargs):
    # [.pyeqn] x = - m * a_2
    result = []
    a_2 = -x / m
    result.append(a_2)
    return result
