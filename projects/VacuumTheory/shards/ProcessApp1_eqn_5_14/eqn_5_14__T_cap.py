from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_5_14__T(self, M: float, P_0: float, W_E: float, **kwargs):
    # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
    result = []
    T = 0.00339889*M*P_0**2/W_E**2
    result.append(T)
    return result
