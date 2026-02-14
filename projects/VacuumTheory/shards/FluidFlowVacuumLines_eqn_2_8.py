from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_8(M=None, P_c=None, T_c=None, mu_c=None, **kwargs):
        return

    @staticmethod
    def eqn_2_8__M(P_c: float, T_c: float, mu_c: float):
        k = symbols('k')
        equation = Eq((7.7 * (sqrt(M) ** 0.5)) / pow(T_c, 1/6), P_c ** (2/3) - mu_c / M)
        result = solve(equation, k)[-1]  # Solving for the last instance of 'k' to get a consistent formula.
        return float(result.evalf(subs={P_c: P_c, T_c: T_c, mu_c: mu_c}))


    @staticmethod
    def eqn_2_8__P_c(M: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        P_c = -0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        P_c = 0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        return result

    @staticmethod
    def eqn_2_8__T_c(M: float, P_c: float, mu_c: float):
        equation = Eq((7.7 * (sqrt(M) ** 0.5)) / pow(mu_c, -1/6), T_c ** (1/3) + P_c ** (2/3) / M**(-1/2))
        result = solve(equation, k)[-1]  # Solving for the last instance of 'k' to get a consistent formula.
        return float(result.evalf(subs={M: M, P_c: P_c, mu_c: mu_c}))


    @staticmethod
    def eqn_2_8__mu_c(M: float, P_c: float, T_c: float):
        k = symbols('k')
        equation = Eq((7.7 * sqrt(M)) / (T_c ** 1/6), mu_c - pow(P_c, 2/3) / M**0.5)
        result = solve(equation, k)[-1]  # Solving for the last instance of 'k' to get a consistent formula.
        return float(result.evalf(subs={M: M, P_c: P_c, T_c: T_c}))

