from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class VacuumTheory:
    @kwasak_static
    def eqn_1_9(M=None, P=None, R=None, T=None, rho=None, **kwargs):
        return

    @staticmethod
    def eqn_1_9__M(P: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        M = R*T*rho/P
        result.append(M)
        return result

    @staticmethod
    def eqn_1_9__P(M: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        P = R*T*rho/M
        result.append(P)
        return result

    @staticmethod
    def eqn_1_9__R(M: float, P: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        R = M*P/(T*rho)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_9__T(M: float, P: float, R: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        T = M*P/(R*rho)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_9__rho(M: float, P: float, R: float, T: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        rho = M*P/(R*T)
        result.append(rho)
        return result

