from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_2__delta(self, lambd: float, psi: float, **kwargs):
    # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
    result = []
    delta = -0.474424998328794*sqrt(lambd/psi)
    result.append(delta)
    delta = 0.474424998328794*sqrt(lambd/psi)
    result.append(delta)
    return result
