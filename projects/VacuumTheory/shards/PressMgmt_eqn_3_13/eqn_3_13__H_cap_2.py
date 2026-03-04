from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_3_13__H_2(self, H_1: float, KAPPA_2: float, P: float, **kwargs):
    # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
    result = []
    H_2 = H_1 + P/KAPPA_2
    result.append(H_2)
    return result
