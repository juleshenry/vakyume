from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_14__T_c(self, T_s: float, **kwargs):
    # [.pyeqn] T_c = T_s + 12
    result = []
    T_c = T_s + 12
    result.append(T_c)
    return result
