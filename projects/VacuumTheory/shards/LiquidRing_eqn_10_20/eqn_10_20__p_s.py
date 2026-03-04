from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_20__p_s(
    self,
    P: float,
    S_0: float,
    S_p: float,
    T_e: float,
    T_i: float,
    p_0: float,
    p_c: float,
    **kwargs,
):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    # Solve for p_s:
    R = (S_0 / (S_p)) ** (1.666666666666667)
    # After clearing **0.6: R = (P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) )
    # (P - p_s) = ((P - p_0)*(460 + T_i) * (P - p_c)) / (R * (P * (460 + T_e)))
    # p_s = P - ((P - p_0)*(460 + T_i) * (P - p_c)) / (R * (P * (460 + T_e)))
    p_s = P - ((P - p_0) * (460 + T_i) * (P - p_c)) / (R * (P * (460 + T_e)))
    return [p_s]
