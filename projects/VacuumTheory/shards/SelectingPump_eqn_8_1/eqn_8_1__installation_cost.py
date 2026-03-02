from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_1__installation_cost(self, NC: float, NS: float, SCON: float, **kwargs):
    # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
    result = []
    installation_cost = 1426.00150101399*SCON**(7/20)*(2.0*NC + NS)
    result.append(installation_cost)
    return result
