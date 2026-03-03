from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_5__P(self, N_i: float, N_nc: float, P_c: float, p_i: float, **kwargs):
    # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
    result = []
    P = P_c + N_nc*p_i/N_i
    result.append(P)
    return result
