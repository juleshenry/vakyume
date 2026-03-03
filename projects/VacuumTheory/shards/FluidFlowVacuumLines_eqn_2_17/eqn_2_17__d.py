from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_17__d(self, L: float, delta_P: float, mu: float, q: float, **kwargs):
    # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
    # Solve for d:
    # Step 1: d ** 4 = delta_P / (0.105 * mu * L * q / 1)
    # Step 2: d = (delta_P / (0.105 * mu * L * q / 1)) ** (1.0 / 4)
    d = (delta_P / (0.105 * mu * L * q / 1)) ** (1.0 / 4)
    return [d]
