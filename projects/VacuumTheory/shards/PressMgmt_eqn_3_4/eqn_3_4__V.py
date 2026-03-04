from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_4__V(self, KAPPA: float, P: float, **kwargs):
    # [.pyeqn] P * V = KAPPA
    result = []
    V = KAPPA/P
    result.append(V)
    return result
