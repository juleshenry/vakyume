from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_5_4__x_i(self, P: float, P_0_i: float, y_i: float, **kwargs):
    # [.pyeqn] y_i * P = x_i * P_0_i
    result = []
    x_i = P * y_i / P_0_i
    result.append(x_i)
    return result
