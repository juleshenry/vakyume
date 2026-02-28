from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_5__t(P_1: float, P_2: float, S_p: float, V: float, **kwargs):
    # [.pyeqn] t = V / S_p * log(P_1 / P_2)
    result = []
    t = V*log(P_1/P_2)/S_p
    result.append(t)
    return result
