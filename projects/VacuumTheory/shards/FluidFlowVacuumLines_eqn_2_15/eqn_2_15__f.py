from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_15__f(self, Re: float, **kwargs):
    # [.pyeqn] f = 0.316 / Re ** (0.25)
    result = []
    f = 0.316/Re**(1/4)
    result.append(f)
    return result
