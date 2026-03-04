from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_5_15__a_M_12(self, M_1: float, M_2: float, P_0_1: float, P_0_2: float, **kwargs):
    # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
    result = []
    a_M_12 = P_0_1*(M_2/M_1)**(2/5)/P_0_2
    result.append(a_M_12)
    return result
