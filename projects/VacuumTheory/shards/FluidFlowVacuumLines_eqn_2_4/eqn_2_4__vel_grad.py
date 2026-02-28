from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_4__vel_grad(_beta: float, mu: float, **kwargs):
    # [.pyeqn] _beta = mu * vel_grad
    result = []
    vel_grad = _beta/mu
    result.append(vel_grad)
    return result
