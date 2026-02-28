from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_17__L(d: float, delta_P: float, mu: float, v: float, **kwargs):
    # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
    result = []
    L = 28.9855072463768*d**2*delta_P/(mu*v)
    result.append(L)
    return result
