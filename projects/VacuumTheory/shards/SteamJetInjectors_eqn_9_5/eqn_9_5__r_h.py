from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_9_5__r_h(self, V: float, t_h: float, w_h: float, **kwargs):
    # [.pyeqn] w_h = r_h * V / t_h
    result = []
    r_h = t_h * w_h / V
    result.append(r_h)
    return result
