from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_12__T_s(self, T_c: float, **kwargs):
    # [.pyeqn] T_c = T_s + 5
    result = []
    T_s = T_c - 5
    result.append(T_s)
    return result
