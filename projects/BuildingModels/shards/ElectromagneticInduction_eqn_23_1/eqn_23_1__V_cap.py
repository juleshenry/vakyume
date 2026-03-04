from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_23_1__V(self, dt: float, **kwargs):
    # [.pyeqn] V = - d%B / dt
    result = []
    V = (Mod(-d, B)) / dt
    result.append(V)
    return result
