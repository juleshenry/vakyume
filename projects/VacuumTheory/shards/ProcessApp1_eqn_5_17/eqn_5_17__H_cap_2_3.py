from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_5_17__H_2_3(
    self, H_2_1: float, H_2_mi: float, x_1: float, x_3: float, **kwargs
):
    # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
    result = []
    H_2_3 = exp((-x_1 * log(H_2_1) + log(H_2_mi)) / x_3)
    result.append(H_2_3)
    return result
