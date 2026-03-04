from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_3__hp(self, installed_costs, **kwargs):
    # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
    result = []
    hp = pow((installed_costs / 38000), 1 / 0.45)
    result.append(hp)
    return [result]
