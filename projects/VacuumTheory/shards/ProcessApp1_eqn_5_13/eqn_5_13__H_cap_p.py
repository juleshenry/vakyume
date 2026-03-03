from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_5_13__H_p(self, HETP: float, N_ES: float, **kwargs):
    # [.pyeqn] H_p = N_ES * HETP
    result = []
    H_p = HETP * N_ES
    result.append(H_p)
    return result
