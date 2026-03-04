from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_2__W_tot(self, F1: float, F2: float, F3: float, x: float, **kwargs):
    # [.pyeqn] W_tot = F1 * x + F2 * x + F3 * x
    result = []
    W_tot = x * (F1 + F2 + F3)
    result.append(W_tot)
    return result
