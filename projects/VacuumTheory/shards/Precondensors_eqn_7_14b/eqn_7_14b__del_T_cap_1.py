from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_7_14b__del_T_1(
    self, A: float, Q_condensor_heat_duty: float, U: float, del_T_2: float, **kwargs
):
    # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
    # del_T_1 appears 2 times — use numerical solver
    from scipy.optimize import brentq
    import numpy as np

    def _res(del_T_1_val):
        try:
            # Force complex evaluation to handle negative bases in fractional powers
            target_var_complex = complex(del_T_1_val, 0)
            val = (Q_condensor_heat_duty / (U * (target_var_complex - del_T_2))) / ln(
                target_var_complex - del_T_2
            ) - A
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
        raise UnsolvedException("No sign change found for del_T_1 in expanded range")
    del_T_1 = brentq(_res, lo, hi)
    return [del_T_1]
