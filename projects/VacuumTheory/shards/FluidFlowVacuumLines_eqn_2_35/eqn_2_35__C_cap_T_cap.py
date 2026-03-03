from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_35__C_T(self, C_L: float, F_p: float, **kwargs):
    # [.pyeqn] C_T = C_L * F_p
    result = []
    C_T = C_L*F_p
    result.append(C_T)
    return result
