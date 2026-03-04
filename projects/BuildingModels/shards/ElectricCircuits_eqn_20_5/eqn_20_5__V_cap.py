from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_20_5__V(self, C: float, IR: float, Q: float, **kwargs):
    # [.pyeqn] V = IR + Q / C
    result = []
    V = IR + Q / C
    result.append(V)
    return result
