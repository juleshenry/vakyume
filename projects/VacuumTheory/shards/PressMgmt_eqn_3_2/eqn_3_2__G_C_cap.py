from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_2__G_C(G: float, H: float, P: float, rho: float, **kwargs):
    # [.pyeqn] P = G / (G_C * rho * H)
    result = []
    G_C = G/(H*P*rho)
    result.append(G_C)
    return result
