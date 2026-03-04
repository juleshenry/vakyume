from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_11_3__F_s(self, t: float, t_c: float, **kwargs):
    # [.pyeqn] t = t_c * F_s
    result = []
    F_s = t/t_c
    result.append(F_s)
    return result
