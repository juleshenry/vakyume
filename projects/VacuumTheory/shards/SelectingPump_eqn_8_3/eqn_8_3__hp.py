from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_3__hp(self, installed_costs: float, **kwargs):
    # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
    # Algebraic: (hp/10)^0.45 = installed_costs/38000
    # hp = 10 * (installed_costs/38000) ^ (1/0.45)
    inner = installed_costs / 38000.0
    hp = 10.0 * inner ** (1.0 / 0.45)
    return [hp]
