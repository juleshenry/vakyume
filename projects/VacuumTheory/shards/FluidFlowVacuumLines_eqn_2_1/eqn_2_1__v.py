from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_1__v(self, D: float, Re: float, mu: float, rho: float, **kwargs):
    # [.pyeqn] Re = rho * D * v / mu
    result = []
    v = Re*mu/(D*rho)
    result.append(v)
    return result
