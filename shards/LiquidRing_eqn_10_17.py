from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_17(
        P: float = None,
        S_0: float = None,
        S_Th: float = None,
        p_0: float = None,
        p_s: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_10_17__P(S_0: float, S_Th: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        try:
            ratio_term = abs(S_Th / S_0) ** (5/3)
            if abs(ratio_term - 1.0) < 1e-12: return []
            P = (p_0 * ratio_term - p_s) / (ratio_term - 1.0)
            result.append(P)
        except:
            pass
        return result

    @staticmethod
    def eqn_10_17__S_0(P: float, S_Th: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        try:
            den = abs((P - p_s) / (P - p_0)) ** 0.6
            if abs(den) < 1e-18: return []
            S_0 = S_Th / den
            result.append(S_0)
        except:
            pass
        return result

    @staticmethod
    def eqn_10_17__S_Th(P: float, S_0: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        try:
            S_Th = S_0 * abs((P - p_s) / (P - p_0)) ** 0.6
            result.append(S_Th)
        except:
            pass
        return result

    @staticmethod
    def eqn_10_17__p_0(P: float, S_0: float, S_Th: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        try:
            ratio_term = abs(S_Th / S_0) ** (5/3)
            if abs(ratio_term) < 1e-18: return []
            p_0 = (P * ratio_term - P + p_s) / ratio_term
            result.append(p_0)
        except:
            pass
        return result

    @staticmethod
    def eqn_10_17__p_s(P: float, S_0: float, S_Th: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        try:
            ratio_term = abs(S_Th / S_0) ** (5/3)
            p_s = P - (P - p_0) * ratio_term
            result.append(p_s)
        except:
            pass
        return result


