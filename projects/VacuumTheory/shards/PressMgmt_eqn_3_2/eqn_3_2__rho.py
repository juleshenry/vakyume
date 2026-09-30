from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_3_2__rho(self, G: float, G_C: float, H: float, P: float, **kwargs):
    # [.pyeqn] P = G / (G_C * rho * H)
    result = []
    rho = G / (G_C * H * P)
    result.append(rho)
    return result
