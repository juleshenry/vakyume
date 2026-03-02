from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_9_5__t_h(self, V: float, r_h: float, w_h: float, **kwargs):
    # [.pyeqn] w_h = r_h * V / t_h
    result = []
    t_h = V*r_h/w_h
    result.append(t_h)
    return result
