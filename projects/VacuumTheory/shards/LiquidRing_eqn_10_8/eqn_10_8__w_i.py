from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_8__w_i(self, bhp: float, c_p: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, **kwargs):
    # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
    result = []
    w_i = 0.02*(-127250.0*bhp + 401.0*c_p*delta_T*f_a*rho)/delta_h_i
    result.append(w_i)
    return result
