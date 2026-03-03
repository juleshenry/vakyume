from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


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
    # k appears transcendentally; use numerical solver
    from scipy.optimize import brentq

    C = w * R * T / (M * 550.0 * 3600.0)
    ratio = float(P_2 / P_1)

    def residual(k_val):
        exponent = (k_val - 1.0) / k_val
        return k_val / (k_val - 1.0) * C * (ratio**exponent - 1.0) - adiabatic_hp

    # Scan a wide range for sign changes
    import numpy as _np

    test_points = _np.concatenate(
        [
            _np.linspace(1.0001, 1.5, 500),
            _np.linspace(1.5, 10.0, 500),
            _np.linspace(10.0, 1000.0, 500),
            _np.linspace(1000.0, 1e6, 500),
        ]
    )
    results = []
    prev_val = None
    for i in range(len(test_points)):
        try:
            cur_val = float(residual(test_points[i]))
        except (ValueError, ZeroDivisionError, OverflowError):
            prev_val = None
            continue
        if prev_val is not None and prev_val * cur_val < 0:
            try:
                k = brentq(residual, float(test_points[i - 1]), float(test_points[i]))
                results.append(k)
            except (ValueError, RuntimeError):
                pass
        prev_val = cur_val

    if results:
        return results
    # No root found — equation may have no solution for k > 1 with these inputs
    raise UnsolvedException("No valid k found in search range")
