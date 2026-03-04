from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_18__P(self, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
    result = []
    P = (S_Th*T_i*p_s + 460*S_Th*p_s - S_p*T_e*p_c - 460*S_p*p_c)/(S_Th*T_i + 460*S_Th - S_p*T_e - 460*S_p)
    result.append(P)
    return result
