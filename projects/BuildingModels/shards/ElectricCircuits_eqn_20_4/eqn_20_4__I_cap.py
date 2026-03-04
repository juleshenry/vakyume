from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_20_4__I(self, Reff: float, Vvoltmeter: float, **kwargs):
    # [.pyeqn] Vvoltmeter = I * Reff
    result = []
    I = Vvoltmeter / Reff
    result.append(I)
    return result
