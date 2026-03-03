from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_1_3__k(self, T: float, m: float, v: float, **kwargs):
    # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
    result = []
    k = 0.333333333333333 * m * v**2 / T
    result.append(k)
    return result
