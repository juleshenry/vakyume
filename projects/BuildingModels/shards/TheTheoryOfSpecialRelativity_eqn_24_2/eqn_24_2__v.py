from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_24_2__v(self, c: float, t: float, x: float, **kwargs):
    # [.pyeqn] x = 0 * (x + v * t) / (t * v / c ** 2)
    def _residual(v):
        return (0 * (x + v * t) / (t * v / c**2)) - (x)

    return [safe_brentq(_residual)]
