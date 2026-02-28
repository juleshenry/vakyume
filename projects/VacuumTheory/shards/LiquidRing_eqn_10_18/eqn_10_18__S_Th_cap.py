from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_18__S_Th(P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
    result = []
    S_Th = S_p*(P*T_e + 460*P - T_e*p_c - 460*p_c)/(P*T_i + 460*P - T_i*p_s - 460*p_s)
    result.append(S_Th)
    return result
