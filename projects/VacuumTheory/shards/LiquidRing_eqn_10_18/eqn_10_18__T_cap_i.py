from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_18__T_i(self, P: float, S_Th: float, S_p: float, T_e: float, p_c: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
    result = []
    T_i = (-460*P*S_Th + P*S_p*T_e + 460*P*S_p + 460*S_Th*p_s - S_p*T_e*p_c - 460*S_p*p_c)/(S_Th*(P - p_s))
    result.append(T_i)
    return result
