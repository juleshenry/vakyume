from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_18__p_c(self, P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
    result = []
    p_c = (-P*S_Th*T_i - 460*P*S_Th + P*S_p*T_e + 460*P*S_p + S_Th*T_i*p_s + 460*S_Th*p_s)/(S_p*(T_e + 460))
    result.append(p_c)
    return result
