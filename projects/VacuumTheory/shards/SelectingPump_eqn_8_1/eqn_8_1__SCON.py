from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_1__SCON(self, NC, NS, installation_cost, **kwargs):
    # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
    result = []
    SCON = pow((installation_cost / (16000 * (NS + 2 * NC))), 1 / 0.35) * 1000
    result.append(SCON)
    return [result]
