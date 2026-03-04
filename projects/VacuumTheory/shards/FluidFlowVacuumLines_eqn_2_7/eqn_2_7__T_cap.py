from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_7__T(self, k: float, m: float, v_a: float, **kwargs):
    # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
    result = []
    T = 0.392699081698724 * m * v_a**2 / k
    result.append(T)
    return result
