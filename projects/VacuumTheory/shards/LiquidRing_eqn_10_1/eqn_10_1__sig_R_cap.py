from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_1__sig_R(self, D_r: float, w: float, **kwargs):
    # [.pyeqn] sig_R = 0.00436 * D_r * w
    result = []
    sig_R = 0.00436 * D_r * w
    result.append(sig_R)
    return result
