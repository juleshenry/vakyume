from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_13__T_c(self, T_s: float, **kwargs):
    # [.pyeqn] T_c = T_s + 25
    result = []
    T_c = T_s + 25
    result.append(T_c)
    return result
