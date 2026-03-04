from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_1__x_0(self, t: float, v_x: float, x: float, **kwargs):
    # [.pyeqn] x = x_0 + v_x * t
    result = []
    x_0 = -t * v_x + x
    result.append(x_0)
    return result
