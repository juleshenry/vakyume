from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_8__rho(self, bhp: float, c_p: float, delta_T: float, delta_h_i: float, f_a: float, w_i: float, **kwargs):
    # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
    result = []
    rho = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*delta_T*f_a)
    result.append(rho)
    return result
