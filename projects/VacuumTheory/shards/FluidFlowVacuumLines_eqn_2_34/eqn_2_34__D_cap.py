from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_2_34__D(
    self, C: float, C_1: float, C_2: float, L: float, P_p: float, mu: float, **kwargs
):
    # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
    # D appears 2 times — use numerical solver
    from scipy.optimize import brentq
    import numpy as np

    def _res(D_val):
        try:
            # Force complex evaluation to handle negative bases in fractional powers
            target_var_complex = complex(D_val, 0)
            val = (
                C_1 * (target_var_complex**4 / (mu * L)) * P_p
                + C_2 * (target_var_complex**3 / L)
                - C
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
        raise UnsolvedException("No sign change found for D in expanded range")
    D = brentq(_res, lo, hi)
    return [D]
