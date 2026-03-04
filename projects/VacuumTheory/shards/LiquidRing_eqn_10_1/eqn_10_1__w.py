from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_1__w(self, D_r: float, sig_R: float, **kwargs):
    # [.pyeqn] sig_R = 0.00436 * D_r * w
    result = []
    w = 229.357798165138*sig_R/D_r
    result.append(w)
    return result
