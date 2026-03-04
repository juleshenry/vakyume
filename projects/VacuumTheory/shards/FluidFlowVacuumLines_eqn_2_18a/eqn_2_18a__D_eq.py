from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_18a__D_eq(self, R_ll: float, **kwargs):
    # [.pyeqn] D_eq = 4 * R_ll
    result = []
    D_eq = 4*R_ll
    result.append(D_eq)
    return result
