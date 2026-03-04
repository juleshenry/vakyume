from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_20_5__IR(self, C: float, Q: float, V: float, **kwargs):
    # [.pyeqn] V = IR + Q / C
    result = []
    IR = V - Q / C
    result.append(IR)
    return result
