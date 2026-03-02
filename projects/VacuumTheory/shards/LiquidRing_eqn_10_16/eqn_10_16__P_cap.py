from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_16__P(self, S_0: float, S_Th: float, p_0: float, **kwargs):
    # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
    result = []
    P = p_0*(S_Th/S_0)**(5/3)/((S_Th/S_0)**1.66666666666667 - 1.0)
    result.append(P)
    P = p_0*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5/((-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
    result.append(P)
    P = p_0*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5/((-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
    result.append(P)
    return result
