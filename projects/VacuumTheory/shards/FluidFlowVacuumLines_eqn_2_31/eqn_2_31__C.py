from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_31__C(self, S_p: float, S_pump_speed: float, **kwargs):
    # [.pyeqn] S_pump_speed = (S_p * C) / (S_p + C)
    result = []
    C = S_p * S_pump_speed / (S_p - S_pump_speed)
    result.append(C)
    return result
