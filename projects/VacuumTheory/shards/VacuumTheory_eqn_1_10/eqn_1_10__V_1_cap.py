from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_10__V_1(self, P_1: float, P_2: float, T_1: float, T_2: float, V_2: float, **kwargs):
    # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
    result = []
    V_1 = P_2*T_1*V_2/(P_1*T_2)
    result.append(V_1)
    return result
