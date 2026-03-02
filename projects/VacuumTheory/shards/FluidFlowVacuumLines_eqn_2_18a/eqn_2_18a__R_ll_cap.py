from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_18a__R_ll(self, D_eq: float, **kwargs):
    # [.pyeqn] D_eq = 4 * R_ll
    result = []
    R_ll = D_eq/4
    result.append(R_ll)
    return result
