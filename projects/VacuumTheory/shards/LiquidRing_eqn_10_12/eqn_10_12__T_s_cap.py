from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_12__T_s(T_c: float, **kwargs):
    # [.pyeqn] T_c = T_s + 5
    result = []
    T_s = T_c - 5
    result.append(T_s)
    return result
