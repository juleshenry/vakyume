from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_17__h_c(R_0: float, R_nc: float, **kwargs):
    # [.pyeqn] R_0 = R_nc + 1 / h_c
    result = []
    h_c = 1/(R_0 - R_nc)
    result.append(h_c)
    return result
