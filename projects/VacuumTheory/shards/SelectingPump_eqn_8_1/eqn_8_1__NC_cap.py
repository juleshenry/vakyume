from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_1__NC(NS: float, SCON: float, installation_cost: float, **kwargs):
    # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
    result = []
    NC = -0.5*NS + 0.000350630766969363*installation_cost/SCON**(7/20)
    result.append(NC)
    return result
