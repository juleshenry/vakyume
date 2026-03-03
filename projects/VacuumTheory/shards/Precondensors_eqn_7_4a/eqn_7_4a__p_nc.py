from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_4a__p_nc(self, P: float, p_c: float, **kwargs):
    # [.pyeqn] p_nc = P - p_c
    result = []
    p_nc = P - p_c
    result.append(p_nc)
    return result
