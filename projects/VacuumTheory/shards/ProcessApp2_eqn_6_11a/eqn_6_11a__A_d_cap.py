from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_11a__A_d(self, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float, **kwargs):
    # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
    result = []
    A_d = delta_h_i*delta_m*m_b/(delta_T*h_d*t_R)
    result.append(A_d)
    return result
