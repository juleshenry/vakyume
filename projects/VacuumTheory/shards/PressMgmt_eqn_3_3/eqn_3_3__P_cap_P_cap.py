from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_3__P_P(self, H_1: float, H_2: float, P: float, **kwargs):
    # [.pyeqn] P_P - P = H_2 - H_1
    result = []
    P_P = -H_1 + H_2 + P
    result.append(P_P)
    return result
