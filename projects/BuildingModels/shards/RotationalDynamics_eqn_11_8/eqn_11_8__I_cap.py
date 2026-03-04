from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_11_8__I(self, i: float, m: float, r: float, **kwargs):
    # [.pyeqn] I = (m * r ** 2) / i
    result = []
    I = m * r**2 / i
    result.append(I)
    return result
