from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_11_6__S_D(
    self,
    P_0_V: float,
    P_D: float,
    P_v_0: float,
    S_B: float,
    p_b: float,
    p_g: float,
    p_v_max: float,
    **kwargs,
):
    # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
    result = []
    S_D = P_D * S_B * (P_0_V - p_b) / (P_D * p_v_max - P_v_0 * p_g - P_v_0 * p_v_max)
    result.append(S_D)
    return result
