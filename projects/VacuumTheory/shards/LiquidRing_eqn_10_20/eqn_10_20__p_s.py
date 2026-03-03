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
    # After raising to 5/3: K = (P-p_0)*(P-p_c) / (P*(P-p_s))
    # p_s is linear: P*(P-p_s) = (P-p_0)*(P-p_c) / K
    # P-p_s = (P-p_0)*(P-p_c) / (K*P)
    # p_s = P - (P-p_0)*(P-p_c) / (K*P)
    K = (S_0 / S_p) ** (5.0 / 3.0) * (460.0 + T_e) / (460.0 + T_i)
    p_s = P - (P - p_0) * (P - p_c) / (K * P)
    return [p_s]
