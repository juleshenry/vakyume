from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_15_8__R(self, R1: float, R2: float, **kwargs):
    # [.pyeqn] R = R1 + R2
    result = []
    R = R1 + R2
    result.append(R)
    return result
