from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_6_10__r_c(self, A: float, dV_dt: float, delta_P: float, mu: float, s: float, tau: float, **kwargs):
    # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
    result = []
    r_c = A*delta_P**(1 - s)/(dV_dt*mu*tau)
    result.append(r_c)
    return result
