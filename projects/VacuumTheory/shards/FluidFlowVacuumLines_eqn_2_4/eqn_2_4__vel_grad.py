from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_4__vel_grad(self, _beta: float, mu: float, **kwargs):
    # [.pyeqn] _beta = mu * vel_grad
    result = []
    vel_grad = _beta / mu
    result.append(vel_grad)
    return result
