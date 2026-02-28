from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_19a__v(R_ll: float, Re: float, mu: float, rho: float, **kwargs):
    # [.pyeqn] Re = 4 * R_ll * rho * v / mu
    result = []
    v = Re*mu/(4*R_ll*rho)
    result.append(v)
    return result
