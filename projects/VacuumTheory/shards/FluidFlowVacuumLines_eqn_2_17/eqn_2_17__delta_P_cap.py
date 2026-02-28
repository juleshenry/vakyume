from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_17__delta_P(L: float, d: float, mu: float, v: float, **kwargs):
    # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
    result = []
    delta_P = 0.0345*L*mu*v/d**2
    result.append(delta_P)
    return result
