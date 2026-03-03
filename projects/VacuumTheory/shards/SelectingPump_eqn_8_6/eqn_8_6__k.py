from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_8_6__k( self, M: float, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, w: float, **kwargs, ):
    # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
    # k appears in the exponent — use numerical solver
    from scipy.optimize import brentq
    def _res(k_val):
        return (k_val / (k_val - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k_val - 1) / k_val) - 1)) - adiabatic_hp
    # Scan for sign change in a wide range
    lo, hi = None, None
    prev = _res(1.01)
    for i in range(1, 10000):
        x = 1.0 + i * 0.01
        try:
            cur = _res(x)
        except Exception:
            continue
        if prev * cur < 0:
            lo, hi = x - 0.01, x
            break
        prev = cur
    if lo is None:
        raise UnsolvedException("No sign change found for k")
    k = brentq(_res, lo, hi)
    return [k]
