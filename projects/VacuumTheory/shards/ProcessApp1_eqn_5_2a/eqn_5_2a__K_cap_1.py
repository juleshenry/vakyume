from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_5_2a__K_1(self, K_2: float, alpha_1_2: float, **kwargs):
    # [.pyeqn] alpha_1_2 = K_1 / K_2
    result = []
    K_1 = K_2*alpha_1_2
    result.append(K_1)
    return result
