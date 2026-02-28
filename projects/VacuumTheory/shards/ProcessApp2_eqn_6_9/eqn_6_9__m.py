from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_9__m(A: float, dV_dt: float, delta_P: float, mu: float, r: float, r_M: float, **kwargs):
    # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
    result = []
    m = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*mu*r)
    result.append(m)
    return result
