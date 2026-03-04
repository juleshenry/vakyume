from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_12_11__L(self, I: float, v: float, **kwargs):
    # [.pyeqn] L = (1 / 2) * I * v ** 2
    result = []
    L = I * v**2 / 2
    result.append(L)
    return result
