from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_6_10__delta_P(self, A: float, dV_dt: float, mu: float, r_c: float, s: float, tau: float, **kwargs):
    # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
    result = []
    delta_P = (dV_dt*mu*r_c*tau/A)**(-1/(s - 1))
    result.append(delta_P)
    return result
