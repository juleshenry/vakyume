from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_3_12__P(self, H_2: float, KAPPA_1: float, **kwargs):
    # [.pyeqn] P = KAPPA_1 * H_2 ** 2
    result = []
    P = H_2**2 * KAPPA_1
    result.append(P)
    return result
