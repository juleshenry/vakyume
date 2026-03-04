from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_7__P_2(self, P_1: float, adiabatic_hp: float, w: float, **kwargs):
    # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
    def _residual(P_2):
        return ((w / 20) * ((P_2 / P_1) ** 0.286 - 1)) - (adiabatic_hp)

    return [safe_brentq(_residual)]
