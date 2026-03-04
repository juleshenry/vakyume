from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_15__P(self, S_Th: float, S_p: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * (P - p_s) / P
    result = []
    P = S_Th*p_s/(S_Th - S_p)
    result.append(P)
    return result
