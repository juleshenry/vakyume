from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_20_4__Reff(self, I: float, Vvoltmeter: float, **kwargs):
    # [.pyeqn] Vvoltmeter = I * Reff
    result = []
    Reff = Vvoltmeter / I
    result.append(Reff)
    return result
