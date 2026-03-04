from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_2__ax(self, t: float, v: float, v_0x: float, **kwargs):
    # [.pyeqn] v = v_0x + ax * t
    result = []
    ax = (v - v_0x) / t
    result.append(ax)
    return result
