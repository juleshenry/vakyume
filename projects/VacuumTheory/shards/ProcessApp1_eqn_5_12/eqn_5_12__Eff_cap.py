from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_12__Eff(N_ES: float, N_t: float, T: float, **kwargs):
    # [.pyeqn] N_t = N_ES / Eff ** T
    result = []
    Eff = (N_ES/N_t)**(1/T)
    result.append(Eff)
    return result
