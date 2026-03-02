from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_11__h_r(self, D: float, L: float, f: float, g_c: float, v: float, **kwargs):
    # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
    x = symbols('x')
    eqn = Eq(f*L*x**2/(2*D*g_c), x)
    try:
        solution = solve(eqn, x)
        h_r = float(solution[0])  # Assuming the first root is the physical one if multiple roots exist.
    except IndexError:
        raise UnsolvedException("The equation has no real solutions.")
    return [h_r]