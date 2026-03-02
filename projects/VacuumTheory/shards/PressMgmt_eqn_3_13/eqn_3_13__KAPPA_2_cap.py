from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_13__KAPPA_2(self, H_1: float, H_2: float, P: float, **kwargs):
    # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
    result = []
    KAPPA_2 = -P/(H_1 - H_2)
    result.append(KAPPA_2)
    return result
