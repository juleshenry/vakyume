from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_11_5__P_D(P_0_v: float, p_g: float, p_v_max: float, **kwargs):
    # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
    result = []
    P_D = P_0_v*(p_g + p_v_max)/p_v_max
    result.append(P_D)
    return result
