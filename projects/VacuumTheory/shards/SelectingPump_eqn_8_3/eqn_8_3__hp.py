from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_3__hp(self, installed_costs: float, **kwargs):
    # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
    # Solve for hp:
    # Step 1: (hp / 10) ** 0.45 = installed_costs / (38000)
    # Step 2: hp / 10 = (installed_costs / (38000)) ** (1.0 / 0.45)
    # Step 3: hp = 10 * (installed_costs / (38000)) ** (1.0 / 0.45)
    hp = 10 * (installed_costs / (38000)) ** (1.0 / 0.45)
    return [hp]
