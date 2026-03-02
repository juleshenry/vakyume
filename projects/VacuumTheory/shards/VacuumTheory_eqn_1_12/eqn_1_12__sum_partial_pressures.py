from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_1_12__sum_partial_pressures(self, Total_P: float, **kwargs):
    # [.pyeqn] Total_P = sum_partial_pressures
    result = []
    sum_partial_pressures = Total_P
    result.append(sum_partial_pressures)
    return result
