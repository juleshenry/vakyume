from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_20__S_0(P: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float, **kwargs):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    result = []
    S_0 = S_p*((P**2*T_i + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P*T_e + 460.0*P - T_e*p_s - 460.0*p_s)))**(3/5)
    result.append(S_0)
    return result
