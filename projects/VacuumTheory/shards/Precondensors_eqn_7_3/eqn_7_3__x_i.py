from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_7_3__x_i(self, P_i_0: float, epsilon_i: float, p_i: float, **kwargs):
    # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
    result = []
    x_i = p_i/(P_i_0*epsilon_i)
    result.append(x_i)
    return result
