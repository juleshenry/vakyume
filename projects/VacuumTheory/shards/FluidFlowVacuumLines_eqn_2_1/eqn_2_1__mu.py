from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_1__mu(self, D: float, Re: float, rho: float, v: float, **kwargs):
    # [.pyeqn] Re = rho * D * v / mu
    result = []
    mu = D*rho*v/Re
    result.append(mu)
    return result
