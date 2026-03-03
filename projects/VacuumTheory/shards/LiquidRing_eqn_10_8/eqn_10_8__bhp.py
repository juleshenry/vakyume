from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_8__bhp(self, c_p: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, w_i: float, **kwargs):
    # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
    result = []
    bhp = 0.00315127701375246*c_p*delta_T*f_a*rho - 0.000392927308447937*delta_h_i*w_i
    result.append(bhp)
    return result
