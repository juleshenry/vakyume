from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_15__Re(self, f: float, **kwargs):
    # [.pyeqn] f = 0.316 / Re ** (0.25)
    result = []
    Re = 0.009971220736/f**4
    result.append(Re)
    return result
