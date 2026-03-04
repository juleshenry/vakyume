from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_6_9__r_M(self, A: float, dV_dt: float, delta_P: float, m: float, mu: float, r: float, **kwargs):
    # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
    result = []
    r_M = A*delta_P/dV_dt - delta_P*m*mu*r/A
    result.append(r_M)
    return result
