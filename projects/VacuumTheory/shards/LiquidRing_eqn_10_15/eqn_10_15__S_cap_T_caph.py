from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_15__S_Th(self, P: float, S_p: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * (P - p_s) / P
    result = []
    S_Th = P * S_p / (P - p_s)
    result.append(S_Th)
    return result
