from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_9__T_c(self, T_s: float, delta_T: float, **kwargs):
    # [.pyeqn] T_c = T_s + delta_T
    result = []
    T_c = T_s + delta_T
    result.append(T_c)
    return result
