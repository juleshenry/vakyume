from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_17a__delta_P(self, L: float, d: float, mu: float, v: float, **kwargs):
    # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
    result = []
    delta_P = 0.0345*L*mu*v/d**2
    result.append(delta_P)
    return result
