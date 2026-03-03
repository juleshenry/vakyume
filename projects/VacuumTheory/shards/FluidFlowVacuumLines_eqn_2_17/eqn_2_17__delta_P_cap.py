from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_17__delta_P(self, L: float, d: float, mu: float, q: float, **kwargs):
    # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
    return [0.105 * mu * L * q / d**4]
