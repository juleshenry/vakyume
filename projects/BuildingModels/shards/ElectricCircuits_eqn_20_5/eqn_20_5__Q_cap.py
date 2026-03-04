from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_20_5__Q(self, C: float, IR: float, V: float, **kwargs):
    # [.pyeqn] V = IR + Q / C
    result = []
    Q = C * (-IR + V)
    result.append(Q)
    return result
