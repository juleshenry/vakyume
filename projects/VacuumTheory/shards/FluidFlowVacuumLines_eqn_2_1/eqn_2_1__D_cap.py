from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_1__D(Re: float, mu: float, rho: float, v: float, **kwargs):
    # [.pyeqn] Re = rho * D * v / mu
    result = []
    D = Re*mu/(rho*v)
    result.append(D)
    return result
