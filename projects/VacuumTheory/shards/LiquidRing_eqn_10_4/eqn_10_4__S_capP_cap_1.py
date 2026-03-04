from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_4__SP_1(
    self, Q_gas: float, SP_2: float, S_p: float, V: float, t: float, **kwargs
):
    # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
    result = []
    SP_1 = Q_gas + (-Q_gas + SP_2) * exp(S_p * t / V)
    result.append(SP_1)
    return result
