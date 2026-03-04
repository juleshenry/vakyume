from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_7_5__P_c(self, N_i: float, N_nc: float, P: float, p_i: float, **kwargs):
    # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
    result = []
    P_c = P - N_nc * p_i / N_i
    result.append(P_c)
    return result
