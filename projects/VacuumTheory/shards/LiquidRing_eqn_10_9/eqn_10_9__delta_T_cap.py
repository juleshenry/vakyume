from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_9__delta_T(self, T_c: float, T_s: float, **kwargs):
    # [.pyeqn] T_c = T_s + delta_T
    result = []
    delta_T = T_c - T_s
    result.append(delta_T)
    return result
