from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_19__S_p(P: float, S_Th: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
    result = []
    S_p = S_Th*((P*T_i + 460.0*P - T_i*p_s - 460.0*p_s)/(P*T_e + 460.0*P - T_e*p_c - 460.0*p_c))**(3/5)
    result.append(S_p)
    return result
