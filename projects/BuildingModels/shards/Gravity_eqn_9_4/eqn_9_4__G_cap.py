from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_9_4__G(self, M: float, U: float, m: float, r: float, **kwargs):
    # [.pyeqn] U = - G * M * m / r
    result = []
    G = -U * r / (M * m)
    result.append(G)
    return result
