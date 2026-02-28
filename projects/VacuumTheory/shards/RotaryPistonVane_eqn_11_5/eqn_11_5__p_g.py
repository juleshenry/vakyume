from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_11_5__p_g(P_0_v: float, P_D: float, p_v_max: float, **kwargs):
    # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
    result = []
    p_g = p_v_max*(-P_0_v + P_D)/P_0_v
    result.append(p_g)
    return result
