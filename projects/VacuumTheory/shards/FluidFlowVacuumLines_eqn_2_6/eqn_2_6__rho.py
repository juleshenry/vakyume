from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_6__rho(self, lambd: float, mu: float, v_a: float, **kwargs):
    # [.pyeqn] mu = 0.35 * rho * lambd * v_a
    result = []
    rho = 2.85714285714286*mu/(lambd*v_a)
    result.append(rho)
    return result
