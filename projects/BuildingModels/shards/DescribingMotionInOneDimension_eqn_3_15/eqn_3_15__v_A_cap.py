from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_15__v_A(self, t: float, v: float, v_B: float, **kwargs):
    # [.pyeqn] v = (v_B + v_A) * t
    result = []
    v_A = -v_B + v / t
    result.append(v_A)
    return result
