from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_15__S_Th(P: float, S_p: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * (P - p_s) / P
    result = []
    S_Th = P*S_p/(P - p_s)
    result.append(S_Th)
    return result
