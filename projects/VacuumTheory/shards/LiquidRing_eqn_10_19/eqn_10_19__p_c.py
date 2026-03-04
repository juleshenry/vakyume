from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_19__p_c(self, P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
    result = []
    p_c = (P*T_e*(S_p/S_Th)**(5/3) - P*T_i + 460.0*P*(S_p/S_Th)**(5/3) - 460.0*P + T_i*p_s + 460.0*p_s)/((S_p/S_Th)**(5/3)*(T_e + 460.0))
    result.append(p_c)
    p_c = (P*T_e*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - P*T_i + 460.0*P*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P + T_i*p_s + 460.0*p_s)/((T_e + 460.0)*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5)
    result.append(p_c)
    p_c = (P*T_e*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - P*T_i + 460.0*P*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P + T_i*p_s + 460.0*p_s)/((T_e + 460.0)*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5)
    result.append(p_c)
    return result
