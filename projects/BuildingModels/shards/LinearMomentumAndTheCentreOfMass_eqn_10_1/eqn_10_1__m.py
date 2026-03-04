from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_1__m(self, p: float, v: float, **kwargs):
    # [.pyeqn] p = m * v
    result = []
    m = p / v
    result.append(m)
    return result
