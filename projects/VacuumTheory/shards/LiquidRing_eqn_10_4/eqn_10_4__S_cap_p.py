from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_4__S_p(
    self, Q_gas: float, SP_1: float, SP_2: float, V: float, t: float, **kwargs
):
    # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
    result = []
    S_p = V * log((Q_gas - SP_1) / (Q_gas - SP_2)) / t
    result.append(S_p)
    return result
