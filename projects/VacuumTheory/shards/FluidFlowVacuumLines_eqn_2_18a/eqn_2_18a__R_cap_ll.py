from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_18a__R_ll(self, D_eq: float, **kwargs):
    # [.pyeqn] D_eq = 4 * R_ll
    result = []
    R_ll = D_eq/4
    result.append(R_ll)
    return result
