from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_5_13__N_ES(self, HETP: float, H_p: float, **kwargs):
    # [.pyeqn] H_p = N_ES * HETP
    result = []
    N_ES = H_p/HETP
    result.append(N_ES)
    return result
