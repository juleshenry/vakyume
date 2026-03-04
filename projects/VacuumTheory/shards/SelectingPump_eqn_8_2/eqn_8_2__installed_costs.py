from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_8_2__installed_costs(self, hp: float, **kwargs):
    # [.pyeqn] installed_costs = 33000 * (hp / 10) ** 0.5
    result = []
    installed_costs = 10435.5162785557*sqrt(hp)
    result.append(installed_costs)
    return result
