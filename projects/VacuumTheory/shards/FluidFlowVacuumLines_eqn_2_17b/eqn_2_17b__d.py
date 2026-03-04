from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_17b__d(self, L: float, delta_P: float, mu: float, q: float, **kwargs):
    # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
    result = []
    d = -0.569242509762222*I*(L*mu*q/delta_P)**(1/4)
    result.append(d)
    d = 0.569242509762222*I*(L*mu*q/delta_P)**(1/4)
    result.append(d)
    d = -0.569242509762222*(L*mu*q/delta_P)**(1/4)
    result.append(d)
    d = 0.569242509762222*(L*mu*q/delta_P)**(1/4)
    result.append(d)
    return result
