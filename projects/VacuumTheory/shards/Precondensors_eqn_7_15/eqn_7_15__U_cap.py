from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_7_15__U(self, sum_R: float, **kwargs):
    # [.pyeqn] 1 / U = sum_R
    result = []
    U = 1/sum_R
    result.append(U)
    return result
