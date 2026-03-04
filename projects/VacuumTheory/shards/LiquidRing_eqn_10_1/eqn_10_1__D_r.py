from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_1__D_r(self, sig_R: float, w: float, **kwargs):
    # [.pyeqn] sig_R = 0.00436 * D_r * w
    result = []
    D_r = 229.357798165138 * sig_R / w
    result.append(D_r)
    return result
