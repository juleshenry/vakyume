from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp2:
    @kwasak_static
    def eqn_6_10(
        A: float = None,
        dV_dt: float = None,
        delta_P: float = None,
        mu: float = None,
        r_c: float = None,
        s: float = None,
        tau: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_6_10__A(
        dV_dt: float, delta_P: float, mu: float, r_c: float, s: float, tau: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        A = dV_dt * delta_P ** (s - 1) * mu * r_c * tau
        result.append(A)
        return result

    @staticmethod
    def eqn_6_10__dV_dt(
        A: float, delta_P: float, mu: float, r_c: float, s: float, tau: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        dV_dt = A * delta_P ** (1 - s) / (mu * r_c * tau)
        result.append(dV_dt)
        return result

    @staticmethod
    def eqn_6_10__delta_P(
        A: float, dV_dt: float, mu: float, r_c: float, s: float, tau: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        delta_P = (dV_dt * mu * r_c * tau / A) ** (-1 / (s - 1))
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_6_10__mu(
        A: float, dV_dt: float, delta_P: float, r_c: float, s: float, tau: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        mu = A * delta_P ** (1 - s) / (dV_dt * r_c * tau)
        result.append(mu)
        return result

    @staticmethod
    def eqn_6_10__r_c(
        A: float, dV_dt: float, delta_P: float, mu: float, s: float, tau: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        r_c = A * delta_P ** (1 - s) / (dV_dt * mu * tau)
        result.append(r_c)
        return result

    @staticmethod
    def eqn_6_10__s(
        A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, tau: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        s = log(A * delta_P / (dV_dt * mu * r_c * tau)) / log(delta_P)
        result.append(s)
        return result

    @staticmethod
    def eqn_6_10__tau(
        A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, s: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        tau = A * delta_P ** (1 - s) / (dV_dt * mu * r_c)
        result.append(tau)
        return result


