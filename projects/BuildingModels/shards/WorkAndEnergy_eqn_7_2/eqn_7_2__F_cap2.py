from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_2__F2(self, F1: float, F3: float, W_tot: float, x: float, **kwargs):
    # [.pyeqn] W_tot = F1 * x + F2 * x + F3 * x
    result = []
    F2 = -F1 - F3 + W_tot / x
    result.append(F2)
    return result
