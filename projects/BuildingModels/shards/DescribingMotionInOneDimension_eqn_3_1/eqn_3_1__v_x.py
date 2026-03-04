from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_1__v_x(self, t: float, x: float, x_0: float, **kwargs):
    # [.pyeqn] x = x_0 + v_x * t
    result = []
    v_x = (x - x_0) / t
    result.append(v_x)
    return result
