from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_1__W(self, F: float, d: float, **kwargs):
    # [.pyeqn] W = F * d
    result = []
    W = F * d
    result.append(W)
    return result
