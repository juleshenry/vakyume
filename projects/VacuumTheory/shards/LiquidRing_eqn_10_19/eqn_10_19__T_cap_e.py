from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_19__T_e(self, P: float, S_Th: float, S_p: float, T_i: float, p_c: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
    # Solve for T_e:
    R = (S_p / (S_Th)) ** (1.666667)
    # (460 + T_e) = ((P - p_s)*(460 + T_i)) / (R * ((P - p_c)))
    # T_e = ((P - p_s)*(460 + T_i)) / (R * ((P - p_c))) - 460
    T_e = ((P - p_s)*(460 + T_i)) / (R * ((P - p_c))) - 460
    return [T_e]
