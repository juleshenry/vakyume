from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_29__C(self, S_1: float, S_2: float, **kwargs):
    # [.pyeqn] S_1 ** -1 = S_2 ** -1 + 1 / C
    result = []
    C = -S_1 * S_2 / (S_1 - S_2)
    result.append(C)
    return result
