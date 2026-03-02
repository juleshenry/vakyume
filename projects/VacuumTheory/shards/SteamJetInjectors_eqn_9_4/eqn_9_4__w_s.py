from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_9_4__w_s(self, AEL: float, SC: float, r: float, **kwargs):
    # [.pyeqn] w_s = AEL * r * SC
    result = []
    w_s = AEL*SC*r
    result.append(w_s)
    return result
