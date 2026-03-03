from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_20__P(
    self,
    S_0: float,
    S_p: float,
    T_e: float,
    T_i: float,
    p_0: float,
    p_c: float,
    p_s: float,
    **kwargs,
):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    # P appears 4 times — use numerical solver
    from scipy.optimize import brentq

    def _res(P_val):
        return (
            S_p
            * (
                (P_val - p_0)
                * (460 + T_i)
                * (P_val - p_c)
                / (P_val * (P_val - p_s) * (460 + T_e))
            )
            ** 0.6
            - S_0
        )

    lo, hi = None, None
    prev = _res(0.01)
    for i in range(1, 100000):
        x = i * 0.01
        try:
            cur = _res(x)
        except Exception:
            continue
        if prev * cur < 0:
            lo, hi = x - 0.01, x
            break
        prev = cur
    if lo is None:
        raise UnsolvedException("No sign change found for P")
    P = brentq(_res, lo, hi)
    return [P]
