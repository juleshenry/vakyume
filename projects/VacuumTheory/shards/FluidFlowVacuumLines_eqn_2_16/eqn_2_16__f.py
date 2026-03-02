from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_16__f(self, Re: float, **kwargs):
    # [.pyeqn] f = 64 / Re
    result = []
    f = 64/Re
    result.append(f)
    return result
