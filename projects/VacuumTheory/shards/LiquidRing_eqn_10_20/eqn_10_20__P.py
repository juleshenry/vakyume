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
    import numpy as np

    def _res(P_val):
        try:
            # Force complex evaluation to handle negative bases in fractional powers
            target_var_complex = complex(P_val, 0)
            val = (
                S_p
                * (
                    (target_var_complex - p_0)
                    * (460 + T_i)
                    * (target_var_complex - p_c)
                    / (target_var_complex * (target_var_complex - p_s) * (460 + T_e))
                )
                ** 0.6
                - S_0
            )
            return val.real if hasattr(val, "real") else val
        except Exception:
            return float("nan")

    lo, hi = None, None
    # Expanded search: log-space from 1e-6 to 1e6 plus some linear steps
    search_points = np.logspace(-6, 6, 500)
    for i in range(len(search_points) - 1):
        p1, p2 = search_points[i], search_points[i + 1]
        r1, r2 = _res(p1), _res(p2)
        if np.isfinite(r1) and np.isfinite(r2) and r1 * r2 <= 0:
            lo, hi = p1, p2
            break
    if lo is None:
        # Fallback to a wider linear search if logspace fails
        for x in np.linspace(0.001, 10000, 1000):
            r = _res(x)
            if np.isfinite(r):
                if lo is None:
                    lo_val, lo = r, x
                if r * lo_val <= 0:
                    hi = x
                    break
    if lo is None or hi is None:
        raise UnsolvedException("No sign change found for P in expanded range")
    P = brentq(_res, lo, hi)
    return [P]
