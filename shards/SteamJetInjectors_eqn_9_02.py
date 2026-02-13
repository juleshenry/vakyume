from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SteamJetInjectors:
    @kwasak_static
    def eqn_9_02(
        P_m: float = None,
        d_n: float = None,
        rho_s: float = None,
        w_s: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_9_02__P_m(d_n: float, rho_s: float, w_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        P_m = 1.334027668054e-6 * w_s**2 / (d_n**4 * rho_s)
        result.append(P_m)
        return result

    @staticmethod
    def eqn_9_02__d_n(P_m: float, rho_s: float, w_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        try:
            d_n = 0.0339853079285911 * sqrt(w_s / (P_m * rho_s) ** 0.5)
            result.append(d_n)
        except:
            pass
        return result

    @staticmethod
    def eqn_9_02__rho_s(P_m: float, d_n: float, w_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        rho_s = 1.334027668054e-6 * w_s**2 / (P_m * d_n**4)
        result.append(rho_s)
        return result

    @staticmethod
    def eqn_9_02__w_s(P_m: float, d_n: float, rho_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        w_s = 865.8 * d_n**2 * sqrt(P_m * rho_s)
        result.append(w_s)
        return result


