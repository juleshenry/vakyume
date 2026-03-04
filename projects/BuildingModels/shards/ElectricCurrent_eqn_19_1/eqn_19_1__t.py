from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_19_1__t(self, I: float, Q: float, **kwargs):
    # [.pyeqn] I = Q / t
    result = []
    t = Q / I
    result.append(t)
    return result
