from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_17__R_0(self, R_nc: float, h_c: float, **kwargs):
    # [.pyeqn] R_0 = R_nc + 1 / h_c
    result = []
    R_0 = R_nc + 1/h_c
    result.append(R_0)
    return result
