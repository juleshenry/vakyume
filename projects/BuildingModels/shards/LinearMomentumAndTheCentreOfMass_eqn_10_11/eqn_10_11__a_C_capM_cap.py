from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_11__a_CM(self, F_ext: float, M: float, **kwargs):
    # [.pyeqn] F_ext = M * a_CM
    result = []
    a_CM = F_ext/M
    result.append(a_CM)
    return result
