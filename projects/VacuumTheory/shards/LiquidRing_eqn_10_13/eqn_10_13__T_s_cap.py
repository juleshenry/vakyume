from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_13__T_s(self, T_c: float, **kwargs):
    # [.pyeqn] T_c = T_s + 25
    result = []
    T_s = T_c - 25
    result.append(T_s)
    return result
