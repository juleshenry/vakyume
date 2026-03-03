from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_3_12__H_2(self, KAPPA_1: float, P: float, **kwargs):
    # [.pyeqn] P = KAPPA_1 * H_2 ** 2
    result = []
    H_2 = -sqrt(P / KAPPA_1)
    result.append(H_2)
    H_2 = sqrt(P / KAPPA_1)
    result.append(H_2)
    return result
