from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_13__H_1(self, H_2: float, KAPPA_2: float, P: float, **kwargs):
    # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
    result = []
    H_1 = H_2 - P/KAPPA_2
    result.append(H_1)
    return result
