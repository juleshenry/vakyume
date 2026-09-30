from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_15__p_s(self, P: float, S_Th: float, S_p: float, **kwargs):
    # [.pyeqn] S_p = S_Th * (P - p_s) / P
    result = []
    p_s = P * (S_Th - S_p) / S_Th
    result.append(p_s)
    return result
