from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_20_5__C(self, IR: float, Q: float, V: float, **kwargs):
    # [.pyeqn] V = IR + Q / C
    result = []
    C = -Q / (IR - V)
    result.append(C)
    return result
