from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_4__t(
    self, Q_gas: float, SP_1: float, SP_2: float, S_p: float, V: float, **kwargs
):
    # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
    result = []
    t = V * log((Q_gas - SP_1) / (Q_gas - SP_2)) / S_p
    result.append(t)
    return result
