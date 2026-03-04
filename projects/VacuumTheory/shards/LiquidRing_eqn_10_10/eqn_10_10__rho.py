from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_10__rho(self, bhp, bhp_0, mu, **kwargs):
    # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
    result = []
    term = log(bhp / bhp_0) - 0.5
    result.append(pow(10, term))
