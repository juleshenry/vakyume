from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_18b__w(R_ll: float, h: float, **kwargs):
    # [.pyeqn] R_ll = w * h / (2 * (w + h))
    result = []
    w = 2*R_ll*h/(-2*R_ll + h)
    result.append(w)
    return result
