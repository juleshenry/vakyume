from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_13__T_s(self, T_c: float, **kwargs):
    # [.pyeqn] T_c = T_s + 25
    result = []
    T_s = T_c - 25
    result.append(T_s)
    return result
