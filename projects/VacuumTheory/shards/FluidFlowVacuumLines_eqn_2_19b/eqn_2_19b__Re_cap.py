from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_19b__Re(h: float, mu: float, rho: float, v: float, w: float, **kwargs):
    # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
    result = []
    Re = 2*h*rho*v*w/(mu*(h + w))
    result.append(Re)
    return result
