from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_19_6__Reff(self, R1: float, R2: float, **kwargs):
    # [.pyeqn] Reff = 1 / (1 / R1 + 1 / R2 + ...)
    def _residual(Reff):
        return (1 / (1 / R1 + 1 / R2 + ...)) - (Reff)

    return [safe_brentq(_residual)]
