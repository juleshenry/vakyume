from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_18b__w(self, R_ll: float, h: float, **kwargs):
    # [.pyeqn] R_ll = w * h / (2 * (w + h))
    result = []
    w = 2 * R_ll * h / (-2 * R_ll + h)
    result.append(w)
    return result
