from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_8_6__k(
    self,
    M: float,
    P_1: float,
    P_2: float,
    R: float,
    T: float,
    adiabatic_hp: float,
    w: float,
    **kwargs,
):
    # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
    def _residual(k):
        return (
            (
                k
                / (k - 1)
                * (w * R * T)
                / (M * 550 * 3600)
                * ((P_2 / P_1) ** ((k - 1) / k) - 1)
            )
        ) - (adiabatic_hp)

    lo = max(0, 1, 1) + 0.01
    hi_candidates = [lo + 0.5, lo + 2, lo + 10, lo * 100, lo + 1000, 1e6]
    from scipy.optimize import brentq as _brentq
    import math as _math

    for hi in hi_candidates:
        try:
            fa = _residual(lo)
            fb = _residual(hi)
            if isinstance(fa, complex):
                fa = fa.real
            if isinstance(fb, complex):
                fb = fb.real
            if _math.isfinite(fa) and _math.isfinite(fb) and fa * fb < 0:

                def _rf(x):
                    v = _residual(x)
                    return v.real if isinstance(v, complex) else float(v)

                return [_brentq(_rf, lo, hi)]
        except Exception:
            continue
    for _lo2, _hi2 in [
        (lo, lo * 1e3),
        (lo, lo * 1e6),
        (lo, lo * 1e8),
        (lo * 0.5 + 0.001, lo * 1e4),
    ]:
        try:
            _fa2 = _residual(_lo2)
            _fb2 = _residual(_hi2)
            if isinstance(_fa2, complex):
                _fa2 = _fa2.real
            if isinstance(_fb2, complex):
                _fb2 = _fb2.real
            if _math.isfinite(_fa2) and _math.isfinite(_fb2) and _fa2 * _fb2 < 0:

                def _rf2(x):
                    v = _residual(x)
                    return v.real if isinstance(v, complex) else float(v)

                return [_brentq(_rf2, _lo2, _hi2)]
        except Exception:
            continue
    return [safe_brentq(_residual)]
