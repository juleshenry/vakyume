from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_5_16__H_i(self, p_i: float, x_i: float, **kwargs):
    # [.pyeqn] p_i = x_i * H_i
    result = []
    H_i = p_i/x_i
    result.append(H_i)
    return result
