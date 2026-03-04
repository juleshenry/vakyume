from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_18b__R_ll(self, h: float, w: float, **kwargs):
    # [.pyeqn] R_ll = w * h / (2 * (w + h))
    result = []
    R_ll = h*w/(2*(h + w))
    result.append(R_ll)
    return result
