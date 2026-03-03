from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_9__mu(self, A: float, dV_dt: float, delta_P: float, m: float, r: float, r_M: float, **kwargs):
    # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
    result = []
    mu = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*r)
    result.append(mu)
    return result
