from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_15__Re(self, f: float, **kwargs):
    # [.pyeqn] f = 0.316 / Re ** (0.25)
    result = []
    Re = 0.009971220736 / f**4
    result.append(Re)
    return result
