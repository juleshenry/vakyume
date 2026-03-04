from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_4_4__x(self, a: float, ax: float, ay: float, y: float, **kwargs):
    # [.pyeqn] a = ax * x + ay * y
    result = []
    x = (a - ay * y) / ax
    result.append(x)
    return result
