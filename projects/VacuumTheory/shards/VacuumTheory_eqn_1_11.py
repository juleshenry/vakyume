from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class VacuumTheory:
    @kwasak_static
    def eqn_1_11(M=None, P=None, T=None, W=None, q=None, **kwargs):
        return

    @staticmethod
    def eqn_1_11__M(P: float, T: float, W: float, q: float, **kwargs):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        M = 6821*T*W/(738*P*q)
        result.append(M)
        return result

    @staticmethod
    def eqn_1_11__P(M: float, T: float, W: float, q: float, **kwargs):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        P = 6821*T*W/(738*M*q)
        result.append(P)
        return result

    @staticmethod
    def eqn_1_11__T(M: float, P: float, W: float, q: float, **kwargs):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        T = 738*M*P*q/(6821*W)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_11__W(M: float, P: float, T: float, q: float, **kwargs):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        W = 738*M*P*q/(6821*T)
        result.append(W)
        return result

    @staticmethod
    def eqn_1_11__q(M: float, P: float, T: float, W: float, **kwargs):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        q = 6821*T*W/(738*M*P)
        result.append(q)
        return result

