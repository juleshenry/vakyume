from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_5__N_i(self, N_nc: float, P: float, P_c: float, p_i: float, **kwargs):
    # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
    result = []
    N_i = N_nc*p_i/(P - P_c)
    result.append(N_i)
    return result
