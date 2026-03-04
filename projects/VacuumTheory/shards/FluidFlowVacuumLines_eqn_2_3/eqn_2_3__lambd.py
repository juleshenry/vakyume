from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_3__lambd(self, D: float, kn: float, **kwargs):
    # [.pyeqn] kn = lambd / D
    result = []
    lambd = D * kn
    result.append(lambd)
    return result
