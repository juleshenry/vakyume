from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_8__c_p(self, bhp: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, w_i: float, **kwargs):
    # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
    result = []
    c_p = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(delta_T*f_a*rho)
    result.append(c_p)
    return result
