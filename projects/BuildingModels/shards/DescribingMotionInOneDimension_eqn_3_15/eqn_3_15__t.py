from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_15__t(self, v: float, v_A: float, v_B: float, **kwargs):
    # [.pyeqn] v = (v_B + v_A) * t
    result = []
    t = v / (v_A + v_B)
    result.append(t)
    return result
