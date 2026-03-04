from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_1__t(self, v_x: float, x: float, x_0: float, **kwargs):
    # [.pyeqn] x = x_0 + v_x * t
    result = []
    t = (x - x_0) / v_x
    result.append(t)
    return result
