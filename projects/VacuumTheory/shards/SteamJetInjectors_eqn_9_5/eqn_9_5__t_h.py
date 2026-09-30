from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_9_5__t_h(self, V: float, r_h: float, w_h: float, **kwargs):
    # [.pyeqn] w_h = r_h * V / t_h
    result = []
    t_h = V * r_h / w_h
    result.append(t_h)
    return result
