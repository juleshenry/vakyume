from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_1__rho(self, D: float, Re: float, mu: float, v: float, **kwargs):
    # [.pyeqn] Re = rho * D * v / mu
    result = []
    rho = Re*mu/(D*v)
    result.append(rho)
    return result
