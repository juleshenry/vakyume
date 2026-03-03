from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_19__S_Th(self, P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
    # [.pyeqn] Rearrange the equation to solve for S_Th based on given inputs and return its real part as a solution
    numerator = P * (T_e - 460) + T_i*(P - p_c)*(460 + T_i) / ((p_s - p_c)*(460 + T_e))**0.6
    denominator = (P - p_c) * (460 + T_e) ** 0.6
    S_Th = numerator / denominator
    return [S_Th]
