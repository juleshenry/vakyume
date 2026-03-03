from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_11a__m_b(self, A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, t_R: float, **kwargs):
    # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
    result = []
    m_b = A_d*delta_T*h_d*t_R/(delta_h_i*delta_m)
    result.append(m_b)
    return result
