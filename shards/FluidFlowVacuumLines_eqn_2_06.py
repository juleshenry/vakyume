from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_06(
        lambd: float = None,
        mu: float = None,
        rho: float = None,
        v_a: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_2_06__lambd(mu: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        lambd = 2.85714285714286 * mu / (rho * v_a)
        result.append(lambd)
        return result

    @staticmethod
    def eqn_2_06__mu(lambd: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        mu = 0.35 * lambd * rho * v_a
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_06__rho(lambd: float, mu: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        rho = 2.85714285714286 * mu / (lambd * v_a)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_06__v_a(lambd: float, mu: float, rho: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        v_a = 2.85714285714286 * mu / (lambd * rho)
        result.append(v_a)
        return result


