from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_31__S_pump_speed(self, C: float, S_p: float, **kwargs):
    # [.pyeqn] S_pump_speed = (S_p * C) / (S_p + C)
    result = []
    S_pump_speed = C*S_p/(C + S_p)
    result.append(S_pump_speed)
    return result
