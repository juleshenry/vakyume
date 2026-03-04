from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_1__rho(self, D: float, Re: float, mu: float, v: float, **kwargs):
    # [.pyeqn] Re = rho * D * v / mu
    result = []
    rho = Re*mu/(D*v)
    result.append(rho)
    return result
