from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_7_5__p_i(self, N_i: float, N_nc: float, P: float, P_c: float, **kwargs):
    # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
    result = []
    p_i = N_i*(P - P_c)/N_nc
    result.append(p_i)
    return result
