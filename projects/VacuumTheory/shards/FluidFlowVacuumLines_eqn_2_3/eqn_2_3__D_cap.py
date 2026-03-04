from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_3__D(self, kn: float, lambd: float, **kwargs):
    # [.pyeqn] kn = lambd / D
    result = []
    D = lambd / kn
    result.append(D)
    return result
