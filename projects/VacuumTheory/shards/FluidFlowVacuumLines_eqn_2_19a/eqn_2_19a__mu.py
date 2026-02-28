from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_19a__mu(R_ll: float, Re: float, rho: float, v: float, **kwargs):
    # [.pyeqn] Re = 4 * R_ll * rho * v / mu
    result = []
    mu = 4*R_ll*rho*v/Re
    result.append(mu)
    return result
