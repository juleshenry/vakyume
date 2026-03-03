from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_1__v(self, m: float, p: float, **kwargs):
    # [.pyeqn] p = m * v
    result = []
    v = p/m
    result.append(v)
    return result
