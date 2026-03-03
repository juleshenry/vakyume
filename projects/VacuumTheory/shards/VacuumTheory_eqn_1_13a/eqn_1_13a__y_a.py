from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_1_13a__y_a(self, n: float, n_a: float, **kwargs):
    # [.pyeqn] y_a = n_a / n
    result = []
    y_a = n_a / n
    result.append(y_a)
    return result
