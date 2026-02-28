from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_11__h_r(D: float, L: float, f: float, g_c: float, v: float, **kwargs):
    # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
    
    x = symbols('x')
    eqn = Eq(f * L * x**2 / (D * 2 * g_c) - x)
    sol = solve(eqn, x)
    h_r_solution = lambdify(v, sol[0])
    
    try:
        result = float(h_r_solution.evalf(subs={'x': v}))
    except IndexError:  # This handles the case where there is no real solution for h_r in terms of x
        raise UnsolvedException("The equation has no real solutions")
    
    return result