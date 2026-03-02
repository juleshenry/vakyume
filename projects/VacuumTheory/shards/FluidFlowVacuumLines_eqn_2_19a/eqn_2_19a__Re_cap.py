from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_19a__Re(self, R_ll: float, mu: float, rho: float, v: float, **kwargs):
    # [.pyeqn] Re = 4 * R_ll * rho * v / mu
    result = []
    Re = 4*R_ll*rho*v/mu
    result.append(Re)
    return result
