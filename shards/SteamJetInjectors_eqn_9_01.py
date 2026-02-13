from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SteamJetInjectors:
    @kwasak_static
    def eqn_9_01(
        A: float = None,
        rho_s: float = None,
        v: float = None,
        w_s: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_9_01__A(rho_s: float, v: float, w_s: float):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        A = w_s / (rho_s * v)
        result.append(A)
        return result

    @staticmethod
    def eqn_9_01__rho_s(A: float, v: float, w_s: float):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        rho_s = w_s / (A * v)
        result.append(rho_s)
        return result

    @staticmethod
    def eqn_9_01__v(A: float, rho_s: float, w_s: float):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        v = w_s / (A * rho_s)
        result.append(v)
        return result

    @staticmethod
    def eqn_9_01__w_s(A: float, rho_s: float, v: float):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        w_s = A * rho_s * v
        result.append(w_s)
        return result


