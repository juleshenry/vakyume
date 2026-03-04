from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_7__P_1(self, P_2, adiabatic_hp, w, **kwargs):
    # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
    result = []
    p_1 = pow((adiabatic_hp * 20 / w + 1), 1 / 0.286)
    result.append(p_1)
    return [result]
