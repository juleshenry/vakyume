from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_6__x_0(self, t: float, v_0: float, x: float, **kwargs):
    # [.pyeqn] x = v_0 * t + x_0
    result = []
    x_0 = -t * v_0 + x
    result.append(x_0)
    return result
