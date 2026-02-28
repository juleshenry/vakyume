from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_1__w(D_r: float, sig_R: float, **kwargs):
    # [.pyeqn] sig_R = 0.00436 * D_r * w
    result = []
    w = 229.357798165138*sig_R/D_r
    result.append(w)
    return result
