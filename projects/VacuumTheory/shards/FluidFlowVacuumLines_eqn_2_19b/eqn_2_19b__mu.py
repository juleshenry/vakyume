from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_19b__mu(self, Re: float, h: float, rho: float, v: float, w: float, **kwargs):
    # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
    result = []
    mu = 2 * h * rho * v * w / (Re * (h + w))
    result.append(mu)
    return result
