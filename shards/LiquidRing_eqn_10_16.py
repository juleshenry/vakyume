from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_16(
        P: float = None,
        S_0: float = None,
        S_Th: float = None,
        p_0: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_10_16__P(S_0: float, S_Th: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        try:
            # S_Th/S_0 = abs(P/(P-p_0))**0.6
            # abs(P/(P-p_0)) = abs(S_Th/S_0)**(5/3)
            ratio_term = abs(S_Th / S_0) ** (5/3)
            if abs(ratio_term - 1.0) < 1e-9:
                return []
            P = p_0 * ratio_term / (ratio_term - 1.0)
            result.append(P)
        except:
            pass
        return result

    @staticmethod
    def eqn_10_16__p_0(P: float, S_0: float, S_Th: float):
        # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        try:
            ratio_term = abs(S_Th / S_0) ** (5/3)
            if abs(ratio_term) < 1e-12:
                return []
            p_0 = P - P / ratio_term
            result.append(p_0)
        except:
            pass
        return result

    @staticmethod
    def eqn_10_16__S_0(P: float, S_Th: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        try:
            if abs(P - p_0) < 1e-12:
                return []
            S_0 = S_Th / (abs(P / (P - p_0)) ** 0.6)
            result.append(S_0)
        except:
            pass
        return result

    @staticmethod
    def eqn_10_16__S_Th(P: float, S_0: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        try:
            if abs(P - p_0) < 1e-12:
                return []
            S_Th = S_0 * (abs(P / (P - p_0)) ** 0.6)
            result.append(S_Th)
        except:
            pass
        return result


