from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_2__x(self, F1: float, F2: float, F3: float, W_tot: float, **kwargs):
    # [.pyeqn] W_tot = F1 * x + F2 * x + F3 * x
    result = []
    x = W_tot / (F1 + F2 + F3)
    result.append(x)
    return result
