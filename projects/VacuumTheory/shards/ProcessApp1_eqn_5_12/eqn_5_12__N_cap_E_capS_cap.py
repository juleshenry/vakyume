from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_12__N_ES(self, Eff: float, N_t: float, T: float, **kwargs):
    # [.pyeqn] N_t = N_ES / Eff ** T
    result = []
    N_ES = Eff**T*N_t
    result.append(N_ES)
    return result
