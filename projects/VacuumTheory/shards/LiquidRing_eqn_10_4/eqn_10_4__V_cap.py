from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_4__V(Q_gas: float, SP_1: float, SP_2: float, S_p: float, t: float, **kwargs):
    # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
    result = []
    V = S_p*t/log((Q_gas - SP_1)/(Q_gas - SP_2))
    result.append(V)
    return result
