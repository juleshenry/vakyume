from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_13__q(self, L: float, d: float, delta_P: float, f: float, rho: float, **kwargs):
    # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
    result = []
    q = -0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
    result.append(q)
    q = 0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
    result.append(q)
    return result
