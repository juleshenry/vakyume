from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_20_4__Vvoltmeter(self, I: float, Reff: float, **kwargs):
    # [.pyeqn] Vvoltmeter = I * Reff
    result = []
    Vvoltmeter = I * Reff
    result.append(Vvoltmeter)
    return result
