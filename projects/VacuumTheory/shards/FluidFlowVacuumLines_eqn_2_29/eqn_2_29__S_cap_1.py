from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_29__S_1(self, C: float, S_2: float, **kwargs):
    # [.pyeqn] S_1 ** -1 = S_2 ** -1 + 1 / C
    result = []
    S_1 = C*S_2/(C + S_2)
    result.append(S_1)
    return result
