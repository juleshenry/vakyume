from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_14__T_s(self, T_c: float, **kwargs):
    # [.pyeqn] T_c = T_s + 12
    result = []
    T_s = T_c - 12
    result.append(T_s)
    return result
