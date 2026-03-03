from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_4__mu(self, _beta: float, vel_grad: float, **kwargs):
    # [.pyeqn] _beta = mu * vel_grad
    result = []
    mu = _beta/vel_grad
    result.append(mu)
    return result
