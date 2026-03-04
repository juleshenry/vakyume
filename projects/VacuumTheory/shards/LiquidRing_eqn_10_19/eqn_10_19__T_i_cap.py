from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_19__T_i(
    self,
    P: float,
    S_Th: float,
    S_p: float,
    T_e: float,
    p_c: float,
    p_s: float,
    **kwargs,
):
    # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
    # Solve for T_i:
    R = (S_p / (S_Th)) ** (1.666666666666667)
    # (460 + T_i) = R * ( (P - p_c)*(460 + T_e) ) / ((P - p_s))
    # T_i = R * ( (P - p_c)*(460 + T_e) ) / ((P - p_s)) - 460
    T_i = R * ((P - p_c) * (460 + T_e)) / ((P - p_s)) - 460
    return [T_i]
