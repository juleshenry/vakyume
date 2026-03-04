from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_1__D(self, Re: float, mu: float, rho: float, v: float, **kwargs):
    # [.pyeqn] Re = rho * D * v / mu
    result = []
    D = Re*mu/(rho*v)
    result.append(D)
    return result
