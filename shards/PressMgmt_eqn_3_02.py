from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class PressMgmt:
    @kwasak_static
    def eqn_3_02(
        G: float = None,
        G_C: float = None,
        H: float = None,
        P: float = None,
        rho: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_3_02__G(G_C: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G = G_C * H * P * rho
        result.append(G)
        return result

    @staticmethod
    def eqn_3_02__G_C(G: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G_C = G / (H * P * rho)
        result.append(G_C)
        return result

    @staticmethod
    def eqn_3_02__H(G: float, G_C: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        H = G / (G_C * P * rho)
        result.append(H)
        return result

    @staticmethod
    def eqn_3_02__P(G: float, G_C: float, H: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        P = G / (G_C * H * rho)
        result.append(P)
        return result

    @staticmethod
    def eqn_3_02__rho(G: float, G_C: float, H: float, P: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        rho = G / (G_C * H * P)
        result.append(rho)
        return result


