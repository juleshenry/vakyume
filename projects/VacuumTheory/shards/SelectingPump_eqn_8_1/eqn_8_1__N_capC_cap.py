from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_1__NC(self, NS: float, SCON: float, installation_cost: float, **kwargs):
    # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
    result = []
    NC = -0.5 * NS + 0.000350630766969363 * installation_cost / SCON ** (7 / 20)
    result.append(NC)
    return result
