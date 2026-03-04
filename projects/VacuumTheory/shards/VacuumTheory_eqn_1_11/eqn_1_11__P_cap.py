from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_1_11__P(self, M: float, T: float, W: float, q: float, **kwargs):
    # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
    result = []
    P = 6821*T*W/(738*M*q)
    result.append(P)
    return result
