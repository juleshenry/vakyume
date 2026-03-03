from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_16__p_0(self, P: float, S_0: float, S_Th: float, **kwargs):
    # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
    result = []
    p_0 = P - P/(S_Th/S_0)**(5/3)
    result.append(p_0)
    p_0 = P - P/(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
    result.append(p_0)
    p_0 = P - P/(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
    result.append(p_0)
    return result
