from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_4_3__v(self, vx: float, vy: float, x: float, y: float, **kwargs):
    # [.pyeqn] v = vx * x + vy * y
    result = []
    v = vx * x + vy * y
    result.append(v)
    return result
