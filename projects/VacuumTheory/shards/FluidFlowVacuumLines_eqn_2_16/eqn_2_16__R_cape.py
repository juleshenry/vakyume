from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_16__Re(self, f: float, **kwargs):
    # [.pyeqn] f = 64 / Re
    result = []
    Re = 64 / f
    result.append(Re)
    return result
