from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_10__r_c(A: float, dV_dt: float, delta_P: float, mu: float, s: float, tau: float, **kwargs):
    # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
    result = []
    r_c = A*delta_P**(1 - s)/(dV_dt*mu*tau)
    result.append(r_c)
    return result
