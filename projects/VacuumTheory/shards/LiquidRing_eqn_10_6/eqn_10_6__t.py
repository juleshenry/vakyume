from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_6__t(self, P_1: float, P_2: float, S_a: float, V: float, **kwargs):
    # [.pyeqn] S_a = V / t * log(P_1 / P_2)
    result = []
    t = V * log(P_1 / P_2) / S_a
    result.append(t)
    return result
