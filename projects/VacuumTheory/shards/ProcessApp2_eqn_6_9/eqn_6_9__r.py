from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_9__r(self, A: float, dV_dt: float, delta_P: float, m: float, mu: float, r_M: float, **kwargs):
    # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
    result = []
    r = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*mu)
    result.append(r)
    return result
