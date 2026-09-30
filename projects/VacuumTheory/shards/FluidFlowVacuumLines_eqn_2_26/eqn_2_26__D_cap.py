from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_2_26__D(
    self,
    L: float,
    P_downstream: float,
    P_p: float,
    P_upstream: float,
    mu: float,
    q: float,
    **kwargs,
):
    # [.pyeqn] q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
    result = []
    D = (
        -14953.4878122122
        * I
        * (
            -L
            * mu
            * q
            / (1.22718463030851e15 * P_downstream - 1.22718463030851e15 * P_upstream)
        )
        ** (1 / 4)
    )
    result.append(D)
    D = (
        14953.4878122122
        * I
        * (
            -L
            * mu
            * q
            / (1.22718463030851e15 * P_downstream - 1.22718463030851e15 * P_upstream)
        )
        ** (1 / 4)
    )
    result.append(D)
    D = -14953.4878122122 * (
        -L
        * mu
        * q
        / (1.22718463030851e15 * P_downstream - 1.22718463030851e15 * P_upstream)
    ) ** (1 / 4)
    result.append(D)
    D = 14953.4878122122 * (
        -L
        * mu
        * q
        / (1.22718463030851e15 * P_downstream - 1.22718463030851e15 * P_upstream)
    ) ** (1 / 4)
    result.append(D)
    return result
