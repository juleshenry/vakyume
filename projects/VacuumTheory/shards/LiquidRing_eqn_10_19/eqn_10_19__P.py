from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_19__P(self, S_Th, S_p, T_e, T_i, p_c, p_s, **kwargs):
    # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
    result = []
    P = (S_p / S_Th * ((460 + T_i) * p_c - (460 + T_e) * p_s) ** (1 / 0.6)) ** (5 / 3)
    if P == 0:
        raise ValueError("P cannot be zero")
    result.append(P)
    return [result]
