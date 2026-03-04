from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_1_10__P_2(self, P_1: float, T_1: float, T_2: float, V_1: float, V_2: float, **kwargs):
    # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
    result = []
    P_2 = P_1*T_2*V_1/(T_1*V_2)
    result.append(P_2)
    return result
