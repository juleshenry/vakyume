from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_1_11__W(self, M: float, P: float, T: float, q: float, **kwargs):
    # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
    result = []
    W = 738*M*P*q/(6821*T)
    result.append(W)
    return result
