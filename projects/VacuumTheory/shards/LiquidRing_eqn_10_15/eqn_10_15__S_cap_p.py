from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_15__S_p(self, P: float, S_Th: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * (P - p_s) / P
    result = []
    S_p = S_Th * (P - p_s) / P
    result.append(S_p)
    return result
