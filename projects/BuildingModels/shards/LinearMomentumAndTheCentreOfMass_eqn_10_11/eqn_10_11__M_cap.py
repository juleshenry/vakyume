from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_11__M(self, F_ext: float, a_CM: float, **kwargs):
    # [.pyeqn] F_ext = M * a_CM
    result = []
    M = F_ext/a_CM
    result.append(M)
    return result
