from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_19b__v(self, Re: float, h: float, mu: float, rho: float, w: float, **kwargs):
    # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
    result = []
    v = Re*mu*(h + w)/(2*h*rho*w)
    result.append(v)
    return result
