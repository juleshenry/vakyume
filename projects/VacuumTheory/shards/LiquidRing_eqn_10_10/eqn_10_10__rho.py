from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_10__rho(self, bhp: float, bhp_0: float, mu: float, **kwargs):
    # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
    def _residual(rho):
        return (bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)) - (bhp)
    return [safe_brentq(_residual)]
