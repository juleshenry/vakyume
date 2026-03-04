from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_9_4__SC(self, AEL: float, r: float, w_s: float, **kwargs):
    # [.pyeqn] w_s = AEL * r * SC
    result = []
    SC = w_s/(AEL*r)
    result.append(SC)
    return result
