from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_19a__rho(self, R_ll: float, Re: float, mu: float, v: float, **kwargs):
    # [.pyeqn] Re = 4 * R_ll * rho * v / mu
    result = []
    rho = Re*mu/(4*R_ll*v)
    result.append(rho)
    return result
