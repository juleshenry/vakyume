from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_18__T_e(P: float, S_Th: float, S_p: float, T_i: float, p_c: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
    result = []
    T_e = (P*S_Th*T_i + 460*P*S_Th - 460*P*S_p - S_Th*T_i*p_s - 460*S_Th*p_s + 460*S_p*p_c)/(S_p*(P - p_c))
    result.append(T_e)
    return result
