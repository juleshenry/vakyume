from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_2_31__S_p(self, C: float, S_pump_speed: float, **kwargs):
    # [.pyeqn] S_pump_speed = (S_p * C) / (S_p + C)
    result = []
    S_p = C*S_pump_speed/(C - S_pump_speed)
    result.append(S_p)
    return result
