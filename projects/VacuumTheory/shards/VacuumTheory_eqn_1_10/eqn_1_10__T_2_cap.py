from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_10__T_2(P_1: float, P_2: float, T_1: float, V_1: float, V_2: float, **kwargs):
    # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
    result = []
    T_2 = P_2*T_1*V_2/(P_1*V_1)
    result.append(T_2)
    return result
