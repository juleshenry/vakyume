from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_20__p_c(
    self,
    P: float,
    S_0: float,
    S_p: float,
    T_e: float,
    T_i: float,
    p_0: float,
    p_s: float,
    **kwargs,
):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    # After raising to 5/3: K = (P-p_0)*(P-p_c) / (P*(P-p_s))
    # p_c is linear: (P-p_c) = K*P*(P-p_s) / (P-p_0)
    # p_c = P - K*P*(P-p_s) / (P-p_0)
    K = (S_0 / S_p) ** (5.0 / 3.0) * (460.0 + T_e) / (460.0 + T_i)
    p_c = P - K * P * (P - p_s) / (P - p_0)
    return [p_c]
