from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_10__rho(self, bhp: float, bhp_0: float, mu: float, **kwargs):
    # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
    # Solve for rho:
    # Step 1: bhp / bhp_0 = 0.5 + 0.0155 *rho ** 0.84* mu ** 0.16
    # Step 2: (bhp / bhp_0 - 0.5) = 0.0155 * rho ** 0.84 * mu ** 0.16
    # Step 3: rho ** 0.84 = ((bhp / bhp_0 - 0.5) / (0.0155 * mu ** 0.16))
    # Step 4: rho = ((bhp / bhp_0 - 0.5) / (0.0155 * mu ** 0.16)) ** (1.0 / 0.84)
    rho = ((bhp / bhp_0 - 0.5) / (0.0155 * mu ** 0.16)) ** (1.0 / 0.84)
    return [rho]
