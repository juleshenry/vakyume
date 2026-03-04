from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_2__G_C(self, G: float, H: float, P: float, rho: float, **kwargs):
    # [.pyeqn] P = G / (G_C * rho * H)
    result = []
    G_C = G / (H * P * rho)
    result.append(G_C)
    return result
