from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_4__hp(self, installed_costs: float, **kwargs):
    # [.pyeqn] installed_costs = 26000 * (hp / 10) ** 0.4
    result = []
    hp = -9.1741667595569e-11*installed_costs**(5/2)
    result.append(hp)
    hp = 9.1741667595569e-11*installed_costs**(5/2)
    result.append(hp)
    return result
