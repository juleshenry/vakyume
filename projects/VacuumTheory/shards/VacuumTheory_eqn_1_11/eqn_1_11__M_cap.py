from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_11__M(P: float, T: float, W: float, q: float, **kwargs):
    # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
    result = []
    M = 6821*T*W/(738*P*q)
    result.append(M)
    return result
