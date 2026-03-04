from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_6__v_0(self, t: float, x: float, x_0: float, **kwargs):
    # [.pyeqn] x = v_0 * t + x_0
    result = []
    v_0 = (x - x_0) / t
    result.append(v_0)
    return result
