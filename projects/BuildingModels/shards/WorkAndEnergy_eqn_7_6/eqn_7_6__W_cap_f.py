from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_6__W_f(self, F_: float, d: float, **kwargs):
    # [.pyeqn] W_f = F_ * d
    result = []
    W_f = F_ * d
    result.append(W_f)
    return result
