from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_9_5__r_h(self, V: float, t_h: float, w_h: float, **kwargs):
    # [.pyeqn] w_h = r_h * V / t_h
    result = []
    r_h = t_h*w_h/V
    result.append(r_h)
    return result
