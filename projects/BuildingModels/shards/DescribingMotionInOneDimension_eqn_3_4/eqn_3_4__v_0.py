from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_4__v_0(self, a: float, t: float, v: float, **kwargs):
    # [.pyeqn] v = v_0 + a * t
    result = []
    v_0 = -a * t + v
    result.append(v_0)
    return result
