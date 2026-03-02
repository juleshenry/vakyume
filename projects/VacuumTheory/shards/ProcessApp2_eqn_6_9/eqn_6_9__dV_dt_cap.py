from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_9__dV_dt(self, A: float, delta_P: float, m: float, mu: float, r: float, r_M: float, **kwargs):
    # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
    result = []
    dV_dt = A**2*delta_P/(A*r_M + delta_P*m*mu*r)
    result.append(dV_dt)
    return result
