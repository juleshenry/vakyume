from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_17__L(self, d: float, delta_P: float, mu: float, q: float, **kwargs):
    # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
    L = (delta_P * d**4) / (0.105 * mu * q)
    return [L]
