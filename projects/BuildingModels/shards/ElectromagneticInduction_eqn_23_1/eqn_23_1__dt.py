from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_23_1__dt(self, V: float, **kwargs):
    # [.pyeqn] V = - d%B / dt
    result = []
    dt = (Mod(-d, B)) / V
    result.append(dt)
    return result
