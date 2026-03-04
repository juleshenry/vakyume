from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_1__d(self, F: float, W: float, **kwargs):
    # [.pyeqn] W = F * d
    result = []
    d = W / F
    result.append(d)
    return result
