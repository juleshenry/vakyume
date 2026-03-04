from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_5_10b__V_1(self, L_0: float, R: float, **kwargs):
    # [.pyeqn] L_0 / V_1 = R / (R + 1)
    result = []
    V_1 = L_0 + L_0/R
    result.append(V_1)
    return result
