from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_7_5__N_nc(self, N_i: float, P: float, P_c: float, p_i: float, **kwargs):
    # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
    result = []
    N_nc = N_i * (P - P_c) / p_i
    result.append(N_nc)
    return result
