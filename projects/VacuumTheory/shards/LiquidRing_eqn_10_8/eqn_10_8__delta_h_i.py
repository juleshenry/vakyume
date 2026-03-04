from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_8__delta_h_i(self, bhp: float, c_p: float, delta_T: float, f_a: float, rho: float, w_i: float, **kwargs):
    # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
    result = []
    delta_h_i = 0.02*(-127250.0*bhp + 401.0*c_p*delta_T*f_a*rho)/w_i
    result.append(delta_h_i)
    return result
