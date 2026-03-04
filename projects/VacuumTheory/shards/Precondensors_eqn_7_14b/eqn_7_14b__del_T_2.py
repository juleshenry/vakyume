from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_14b__del_T_2( self, A: float, Q_condensor_heat_duty: float, U: float, del_T_1: float, **kwargs ):
    # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / log(del_T_1 - del_T_2)
    # del_T_2 appears 2 times — use numerical solver
    from scipy.optimize import brentq
    import numpy as np
    import math as _math
    def _res(del_T_2_val):
        try:
            log = _math.log
            sqrt = _math.sqrt
            exp = _math.exp
            val = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2_val))) / log(del_T_1 - del_T_2_val) - A
            return float(val)
        except Exception:
            return float('nan')
    lo, hi = None, None
    # Build search points: logspace + offsets from parameter values
    pos_pts = list(np.logspace(-6, 7, 600))
    neg_pts = [-x for x in reversed(pos_pts)]
    # Add parameter-based offsets to catch roots near param values
    param_vals = [A, Q_condensor_heat_duty, U, del_T_1]
    extra_pts = []
    for pv in param_vals:
        for off in [-100, -10, -1, -0.1, -0.01, 0.01, 0.1, 0.5, 1, 1.01, 2, 10, 100, 1000, 10000]:
            extra_pts.append(pv + off)
    search_points = sorted(set(neg_pts + [0.0] + pos_pts + extra_pts))
    best_lo, best_hi, best_width = None, None, float('inf')
    prev_r = None
    for pt in search_points:
        r = _res(pt)
        if not np.isfinite(r):
            prev_r = None
            continue
        if prev_r is not None and prev_r * r <= 0:
            width = abs(pt - prev_pt)
            if width < best_width:
                best_lo, best_hi, best_width = prev_pt, pt, width
        prev_r, prev_pt = r, pt
    if best_lo is not None:
        lo, hi = best_lo, best_hi
    if lo is None:
        raise UnsolvedException("No sign change found for del_T_2 in expanded range")
    mid = (lo + hi) / 2.0
    mid_r = _res(mid)
    if not np.isfinite(mid_r):
        # Bracket spans a discontinuity; search for a valid sub-bracket
        found_sub = False
        for sub_pts in [np.linspace(lo, mid, 50), np.linspace(mid, hi, 50)]:
            sub_prev_r, sub_prev_pt = None, None
            for sp in sub_pts:
                sr = _res(sp)
                if not np.isfinite(sr):
                    sub_prev_r = None
                    continue
                if sub_prev_r is not None and sub_prev_r * sr <= 0:
                    lo, hi = sub_prev_pt, sp
                    found_sub = True
                    break
                sub_prev_r, sub_prev_pt = sr, sp
            if found_sub:
                break
        if not found_sub:
            raise UnsolvedException("Bracket spans discontinuity for del_T_2")
    del_T_2 = brentq(_res, lo, hi)
    return [del_T_2]
