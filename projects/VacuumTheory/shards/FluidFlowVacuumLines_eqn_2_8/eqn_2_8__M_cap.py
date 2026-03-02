from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_8__M(self, P_c: float, T_c: float, mu_c: float, **kwargs):
    # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
    result = []
    M = 0.0168662506324844*T_c**(1/3)*mu_c**2/P_c**(4/3)
    result.append(M)
    return result
