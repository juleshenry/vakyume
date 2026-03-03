from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_9__T_s(self, T_c: float, delta_T: float, **kwargs):
    # [.pyeqn] T_c = T_s + delta_T
    result = []
    T_s = T_c - delta_T
    result.append(T_s)
    return result
