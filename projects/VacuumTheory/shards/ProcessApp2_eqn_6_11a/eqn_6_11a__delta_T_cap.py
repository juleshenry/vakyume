from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_11a__delta_T(self, A_d: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float, **kwargs):
    # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
    result = []
    delta_T = delta_h_i*delta_m*m_b/(A_d*h_d*t_R)
    result.append(delta_T)
    return result
