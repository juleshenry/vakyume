from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_2__G(self, G_C: float, H: float, P: float, rho: float, **kwargs):
    # [.pyeqn] P = G / (G_C * rho * H)
    result = []
    G = G_C*H*P*rho
    result.append(G)
    return result
