from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_4__Q_gas(
    self, SP_1: float, SP_2: float, S_p: float, V: float, t: float, **kwargs
):
    # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
    result = []
    Q_gas = -(SP_1 - SP_2 * exp(S_p * t / V)) / (exp(S_p * t / V) - 1)
    result.append(Q_gas)
    return result
