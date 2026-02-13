from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp2:
    @kwasak_static
    def eqn_6_09(
        A: float = None,
        dV_dt: float = None,
        delta_P: float = None,
        m: float = None,
        mu: float = None,
        r: float = None,
        r_M: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_6_09__A(
        dV_dt: float, delta_P: float, m: float, mu: float, r: float, r_M: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        try:
            disc = dV_dt**2 * r_M**2 + 4 * delta_P**2 * dV_dt * mu * m * r
            if disc < 0: return []
            A1 = (dV_dt * r_M - sqrt(disc)) / (2 * delta_P)
            A2 = (dV_dt * r_M + sqrt(disc)) / (2 * delta_P)
            # Area is usually positive
            if A1 > 0: result.append(A1)
            if A2 > 0: result.append(A2)
            if not result:
                result = [A1, A2]
        except:
            pass
        return result

    @staticmethod
    def eqn_6_09__dV_dt(
        A: float, delta_P: float, m: float, mu: float, r: float, r_M: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        try:
            den = (A * r_M + delta_P * m * mu * r)
            if abs(den) < 1e-12: return []
            dV_dt = A**2 * delta_P / den
            result.append(dV_dt)
        except:
            pass
        return result

    @staticmethod
    def eqn_6_09__delta_P(
        A: float, dV_dt: float, m: float, mu: float, r: float, r_M: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        try:
            den = (A**2 - dV_dt * m * mu * r)
            if abs(den) < 1e-12: return []
            delta_P = A * dV_dt * r_M / den
            result.append(delta_P)
        except:
            pass
        return result

    @staticmethod
    def eqn_6_09__m(
        A: float, dV_dt: float, delta_P: float, mu: float, r: float, r_M: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        try:
            den = (dV_dt * delta_P * mu * r)
            if abs(den) < 1e-12: return []
            m = A * (A * delta_P - dV_dt * r_M) / den
            result.append(m)
        except:
            pass
        return result

    @staticmethod
    def eqn_6_09__mu(
        A: float, dV_dt: float, delta_P: float, m: float, r: float, r_M: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        try:
            den = (dV_dt * delta_P * m * r)
            if abs(den) < 1e-12: return []
            mu = A * (A * delta_P - dV_dt * r_M) / den
            result.append(mu)
        except:
            pass
        return result

    @staticmethod
    def eqn_6_09__r(
        A: float, dV_dt: float, delta_P: float, m: float, mu: float, r_M: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        try:
            den = (dV_dt * delta_P * m * mu)
            if abs(den) < 1e-12: return []
            r = A * (A * delta_P - dV_dt * r_M) / den
            result.append(r)
        except:
            pass
        return result

    @staticmethod
    def eqn_6_09__r_M(
        A: float, dV_dt: float, delta_P: float, m: float, mu: float, r: float
    ):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        try:
            if abs(dV_dt) < 1e-12 or abs(A) < 1e-12: return []
            r_M = (A * delta_P / dV_dt) - (delta_P * m * mu * r / A)
            result.append(r_M)
        except:
            pass
        return result


