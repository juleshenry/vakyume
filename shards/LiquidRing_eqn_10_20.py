from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_20(
        P: float = None,
        S_0: float = None,
        S_p: float = None,
        T_e: float = None,
        T_i: float = None,
        p_0: float = None,
        p_c: float = None,
        p_s: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_10_20__P(
        S_0: float,
        S_p: float,
        T_e: float,
        T_i: float,
        p_0: float,
        p_c: float,
        p_s: float,
    ):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        def equation(P_val):
            try:
                den = (P_val * (P_val - p_s) * (460 + T_e))
                if abs(den) < 1e-12: return 1e12
                ratio = ((P_val - p_0) * (460 + T_i) * (P_val - p_c)) / den
                return S_0 - S_p * (abs(ratio) ** 0.6)
            except:
                return 1e12

        try:
            from scipy.optimize import fsolve
            # Try multiple guesses to cover different regions
            for guess in [max(p_0, p_c, p_s) + 1.0, 100.0, 1.0, 0.1, 1000.0, -100.0]:
                result, info, ier, msg = fsolve(equation, guess, full_output=True)
                if ier == 1:
                    if abs(equation(result[0])) < 1e-6:
                        return [float(result[0])]
            return []
        except:
            return []

    @staticmethod
    def eqn_10_20__S_0(
        P: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float
    ):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        try:
            den = (P * (P - p_s) * (460 + T_e))
            if abs(den) < 1e-12: return []
            base = ((P - p_0) * (460 + T_i) * (P - p_c)) / den
            S_0 = S_p * (abs(base) ** 0.6)
            result.append(S_0)
        except:
            pass
        return result

    @staticmethod
    def eqn_10_20__S_p(
        P: float, S_0: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float
    ):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        try:
            den = (P * (P - p_s) * (460 + T_e))
            if abs(den) < 1e-12: return []
            base = ((P - p_0) * (460 + T_i) * (P - p_c)) / den
            # S_0 = S_p * abs(base)**0.6 => S_p = S_0 / abs(base)**0.6
            denominator = abs(base) ** 0.6
            if abs(denominator) < 1e-18: return []
            S_p = S_0 / denominator
            result.append(S_p)
        except:
            pass
        return result

    @staticmethod
    def eqn_10_20__T_e(
        P: float, S_0: float, S_p: float, T_i: float, p_0: float, p_c: float, p_s: float
    ):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        try:
            # (S_0/S_p)**(5/3) = ((P-p_0)*(460+T_i)*(P-p_c)) / (P*(P-p_s)*(460+T_e))
            # Let K = (S_0/S_p)**(5/3)
            # 460 + T_e = ((P-p_0)*(460+T_i)*(P-p_c)) / (K * P * (P-p_s))
            K = abs(S_0 / S_p) ** (5/3)
            if abs(K) < 1e-18: return []
            num = (P - p_0) * (460 + T_i) * (P - p_c)
            den = K * P * (P - p_s)
            if abs(den) < 1e-18: return []
            T_e = (num / den) - 460
            return [T_e]
        except:
            return []

    @staticmethod
    def eqn_10_20__T_i(
        P: float, S_0: float, S_p: float, T_e: float, p_0: float, p_c: float, p_s: float
    ):
        try:
            K = abs(S_0 / S_p) ** (5/3)
            num = K * P * (P - p_s) * (460 + T_e)
            den = (P - p_0) * (P - p_c)
            if abs(den) < 1e-18: return []
            T_i = (num / den) - 460
            return [T_i]
        except:
            return []

    @staticmethod
    def eqn_10_20__p_c(
        P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_s: float
    ):
        try:
            K = abs(S_0 / S_p) ** (5/3)
            # K = ((P-p_0)*(460+T_i)*(P-p_c)) / (P*(P-p_s)*(460+T_e))
            # P - p_c = (K * P * (P-p_s) * (460+T_e)) / ((P-p_0)*(460+T_i))
            den = (P - p_0) * (460 + T_i)
            if abs(den) < 1e-18: return []
            p_c = P - (K * P * (P - p_s) * (460 + T_e)) / den
            return [p_c]
        except:
            return []

    @staticmethod
    def eqn_10_20__p_s(
        P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float
    ):
        try:
            K = abs(S_0 / S_p) ** (5/3)
            # K = ((P-p_0)*(460+T_i)*(P-p_c)) / (P*(P-p_s)*(460+T_e))
            # P - p_s = ((P-p_0)*(460+T_i)*(P-p_c)) / (K * P * (460+T_e))
            den = K * P * (460 + T_e)
            if abs(den) < 1e-18: return []
            p_s = P - ((P - p_0) * (460 + T_i) * (P - p_c)) / den
            return [p_s]
        except:
            return []


