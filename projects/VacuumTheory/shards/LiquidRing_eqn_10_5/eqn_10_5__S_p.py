from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_5__S_p(self, P_1: float, P_2: float, V: float, t: float, **kwargs):
    # [.pyeqn] t = V / S_p * log(P_1 / P_2)
    result = []
    S_p = V * log(P_1 / P_2) / t
    result.append(S_p)
    return result
