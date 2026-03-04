from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_12_11__I(self, L: float, v: float, **kwargs):
    # [.pyeqn] L = (1 / 2) * I * v ** 2
    result = []
    I = 2 * L / v**2
    result.append(I)
    return result
