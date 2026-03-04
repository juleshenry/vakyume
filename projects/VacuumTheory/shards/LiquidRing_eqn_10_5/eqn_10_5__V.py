from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_5__V(self, P_1: float, P_2: float, S_p: float, t: float, **kwargs):
    # [.pyeqn] t = V / S_p * log(P_1 / P_2)
    result = []
    V = S_p*t/log(P_1/P_2)
    result.append(V)
    return result
