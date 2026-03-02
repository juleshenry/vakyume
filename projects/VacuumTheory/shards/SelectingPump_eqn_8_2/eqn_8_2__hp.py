from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_2__hp(self, installed_costs: float, **kwargs):
    # [.pyeqn] installed_costs = 33000 * (hp / 10) ** 0.5
    result = []
    hp = 9.18273645546364e-9*installed_costs**2
    result.append(hp)
    return result
