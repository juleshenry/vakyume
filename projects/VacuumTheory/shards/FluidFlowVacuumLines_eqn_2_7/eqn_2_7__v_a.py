from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_7__v_a(self, T: float, k: float, m: float, **kwargs):
    # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
    result = []
    v_a = 1.59576912160573 * sqrt(T * k / m)
    result.append(v_a)
    return result
