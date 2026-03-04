from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_4__t(self, a: float, v: float, v_0: float, **kwargs):
    # [.pyeqn] v = v_0 + a * t
    result = []
    t = (v - v_0) / a
    result.append(t)
    return result
