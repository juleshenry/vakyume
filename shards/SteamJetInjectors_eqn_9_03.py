from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SteamJetInjectors:
    @kwasak_static
    def eqn_9_03(
        P_s: float = None,
        V: float = None,
        t_e: float = None,
        w_j: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_9_03__P_s(V: float, t_e: float, w_j: float):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        P_s = 33.3333333333333 * (23.0 * V - 10.0 * t_e * w_j) / V
        result.append(P_s)
        return result

    @staticmethod
    def eqn_9_03__V(P_s: float, t_e: float, w_j: float):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        V = -1000.0 * t_e * w_j / (3.0 * P_s - 2300.0)
        result.append(V)
        return result

    @staticmethod
    def eqn_9_03__t_e(P_s: float, V: float, w_j: float):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        t_e = 0.001 * V * (2300.0 - 3.0 * P_s) / w_j
        result.append(t_e)
        return result

    @staticmethod
    def eqn_9_03__w_j(P_s: float, V: float, t_e: float):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        w_j = 0.001 * V * (2300.0 - 3.0 * P_s) / t_e
        result.append(w_j)
        return result


