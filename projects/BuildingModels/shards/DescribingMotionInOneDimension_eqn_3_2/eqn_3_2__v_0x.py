from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_2__v_0x(self, ax: float, t: float, v: float, **kwargs):
    # [.pyeqn] v = v_0x + ax * t
    result = []
    v_0x = -ax * t + v
    result.append(v_0x)
    return result
