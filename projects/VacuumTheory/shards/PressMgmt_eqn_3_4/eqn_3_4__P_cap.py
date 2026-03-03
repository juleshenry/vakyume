from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_4__P(self, KAPPA: float, V: float, **kwargs):
    # [.pyeqn] P * V = KAPPA
    result = []
    P = KAPPA/V
    result.append(P)
    return result
