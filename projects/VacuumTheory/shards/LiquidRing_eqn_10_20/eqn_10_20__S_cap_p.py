from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_20__S_p(self, P: float, S_0: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float, **kwargs):
    # [.pyeqn] Solve for S_p by rearranging the equation
    term1 = (P - p_0)**(4/3) * ((T_e + 460) / (T_i + T_e))**(2/3)
    # Rearrange the equation to solve for S_p: S_p = S_0 / term1
    S_p = S_0 / term1
    return [S_p]
