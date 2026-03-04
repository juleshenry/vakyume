from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_11_5__p_v_max(self, P_0_v: float, P_D: float, p_g: float, **kwargs):
    # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
    result = []
    p_v_max = -P_0_v * p_g / (P_0_v - P_D)
    result.append(p_v_max)
    return result
