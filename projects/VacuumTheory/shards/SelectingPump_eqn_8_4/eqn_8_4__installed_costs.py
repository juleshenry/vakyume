from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_4__installed_costs(hp: float, **kwargs):
    # [.pyeqn] installed_costs = 26000 * (hp / 10) ** 0.4
    result = []
    installed_costs = 10350.7864343909*hp**(2/5)
    result.append(installed_costs)
    return result
