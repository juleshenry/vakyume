from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_2__lambd(self, delta: float, psi: float, **kwargs):
    # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
    result = []
    lambd = 4.44288293815837*delta**2*psi
    result.append(lambd)
    return result
