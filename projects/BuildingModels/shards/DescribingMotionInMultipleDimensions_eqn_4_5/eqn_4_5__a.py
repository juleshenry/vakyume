from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_4_5__a(self, t: float, v: float, v_0: float, **kwargs):
    # [.pyeqn] v = v_0 + a * t
    result = []
    a = (v - v_0) / t
    result.append(a)
    return result
