from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_21(P=None, P_d=None, P_prime=None, **kwargs):
        return

    @staticmethod
    def eqn_10_21__P(P_d: float, P_prime: float, **kwargs):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P = P_d*P_prime/760
        result.append(P)
        return result

    @staticmethod
    def eqn_10_21__P_d(P: float, P_prime: float, **kwargs):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P_d = 760*P/P_prime
        result.append(P_d)
        return result

    @staticmethod
    def eqn_10_21__P_prime(P: float, P_d: float, **kwargs):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P_prime = 760*P/P_d
        result.append(P_prime)
        return result

