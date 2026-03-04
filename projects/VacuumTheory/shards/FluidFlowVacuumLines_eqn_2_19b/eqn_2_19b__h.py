from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_19b__h(self, Re: float, mu: float, rho: float, v: float, w: float, **kwargs):
    # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
    result = []
    h = Re*mu*w/(-Re*mu + 2*rho*v*w)
    result.append(h)
    return result
