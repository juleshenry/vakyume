from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_1__SCON(self, NC, NS, installation_cost, **kwargs):
    # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
    def _residual(SCON):
        return (16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35) - (installation_cost)

    return [safe_brentq(_residual)]
