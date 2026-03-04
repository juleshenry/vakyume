from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_7_2__P_i_0(self, p_i: float, x_i: float, **kwargs):
    # [.pyeqn] p_i = x_i * P_i_0
    result = []
    P_i_0 = p_i/x_i
    result.append(P_i_0)
    return result
