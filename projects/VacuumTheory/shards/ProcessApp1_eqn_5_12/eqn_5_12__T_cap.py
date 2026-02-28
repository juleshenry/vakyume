from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_12__T(Eff: float, N_ES: float, N_t: float, **kwargs):
    # [.pyeqn] N_t = N_ES / Eff ** T
    result = []
    T = log(N_ES/N_t)/log(Eff)
    result.append(T)
    return result
