from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_4__SP_2(Q_gas: float, SP_1: float, S_p: float, V: float, t: float, **kwargs):
    # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
    result = []
    SP_2 = (Q_gas*exp(S_p*t/V) - Q_gas + SP_1)*exp(-S_p*t/V)
    result.append(SP_2)
    return result
