from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_2__psi(self, delta: float, lambd: float, **kwargs):
    # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
    result = []
    psi = 0.225079079039277*lambd/delta**2
    result.append(psi)
    return result
