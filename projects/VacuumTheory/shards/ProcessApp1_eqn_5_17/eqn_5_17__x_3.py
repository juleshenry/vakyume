from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_17__x_3(self, H_2_1: float, H_2_3: float, H_2_mi: float, x_1: float, **kwargs):
    # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
    result = []
    x_3 = (-x_1*log(H_2_1) + log(H_2_mi))/log(H_2_3)
    result.append(x_3)
    return result
