from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_17__R_nc(self, R_0: float, h_c: float, **kwargs):
    # [.pyeqn] R_0 = R_nc + 1 / h_c
    result = []
    R_nc = R_0 - 1 / h_c
    result.append(R_nc)
    return result
