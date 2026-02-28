from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_3_2__rho(G: float, G_C: float, H: float, P: float, **kwargs):
    # [.pyeqn] P = G / (G_C * rho * H)
    result = []
    rho = G/(G_C*H*P)
    result.append(rho)
    return result
