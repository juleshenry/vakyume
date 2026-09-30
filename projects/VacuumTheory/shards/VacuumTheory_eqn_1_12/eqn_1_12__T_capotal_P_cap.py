from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_1_12__Total_P(self, sum_partial_pressures: float, **kwargs):
    # [.pyeqn] Total_P = sum_partial_pressures
    result = []
    Total_P = sum_partial_pressures
    result.append(Total_P)
    return result
