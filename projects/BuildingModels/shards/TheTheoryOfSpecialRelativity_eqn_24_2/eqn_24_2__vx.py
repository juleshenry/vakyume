from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_2__vx(self, c: float, t: float, v: float, x: float, **kwargs):
    # [.pyeqn] x = 0 * (x + v * t) + vx / c ** 2
    result = []
    vx = c**2 * x
    result.append(vx)
    return result
