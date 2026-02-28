from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_1__sig_R(D_r: float, w: float, **kwargs):
    # [.pyeqn] sig_R = 0.00436 * D_r * w
    result = []
    sig_R = 0.00436*D_r*w
    result.append(sig_R)
    return result
