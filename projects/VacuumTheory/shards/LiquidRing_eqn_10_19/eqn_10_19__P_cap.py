from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_19__P(
    self,
    S_Th: float,
    S_p: float,
    T_e: float,
    T_i: float,
    p_c: float,
    p_s: float,
    **kwargs,
):
    # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
    # Raise both sides to 5/3: R = (S_p/S_Th)^(5/3) * (460+T_e)/(460+T_i)
    # Then R = (P - p_s)/(P - p_c)  =>  R*(P - p_c) = P - p_s  =>  P*(R-1) = R*p_c - p_s
    # P = (R*p_c - p_s) / (R - 1)
    R = (S_p / S_Th) ** (5.0 / 3.0) * (460.0 + T_e) / (460.0 + T_i)
    P = (R * p_c - p_s) / (R - 1.0)
    return [P]
