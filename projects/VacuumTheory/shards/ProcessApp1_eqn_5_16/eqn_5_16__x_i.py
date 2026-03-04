from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_5_16__x_i(self, H_i: float, p_i: float, **kwargs):
    # [.pyeqn] p_i = x_i * H_i
    result = []
    x_i = p_i / H_i
    result.append(x_i)
    return result
