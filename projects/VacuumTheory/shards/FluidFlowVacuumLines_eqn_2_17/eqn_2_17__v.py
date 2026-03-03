from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_17__v(self, L: float, d: float, delta_P: float, mu: float, **kwargs):
    # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
    return (delta_P * d**2) / (mu * L)
