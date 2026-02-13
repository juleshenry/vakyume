from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_19(
        P: float = None,
        S_Th: float = None,
        S_p: float = None,
        T_e: float = None,
        T_i: float = None,
        p_c: float = None,
        p_s: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_10_19__P(
        S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float
    ):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        try:
            power_term = (S_p / S_Th) ** (5/3)
            # power_term = ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))
            # power_term * (P - p_c) * (460 + T_e) = (P - p_s) * (460 + T_i)
            # P * power_term * (460 + T_e) - p_c * power_term * (460 + T_e) = P * (460 + T_i) - p_s * (460 + T_i)
            # P * (power_term * (460 + T_e) - (460 + T_i)) = p_c * power_term * (460 + T_e) - p_s * (460 + T_i)
            P = (p_c * power_term * (460 + T_e) - p_s * (460 + T_i)) / (power_term * (460 + T_e) - (460 + T_i))
            return [P]
        except:
            return []

    @staticmethod
    def eqn_10_19__S_Th(
        P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float
    ):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        try:
            base = ((P - p_s) * (460 + T_i)) / ((P - p_c) * (460 + T_e))
            S_Th = S_p / (abs(base) ** 0.6)
            return [S_Th]
        except:
            return []

    @staticmethod
    def eqn_10_19__S_p(
        P: float, S_Th: float, T_e: float, T_i: float, p_c: float, p_s: float
    ):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        try:
            base = ((P - p_s) * (460 + T_i)) / ((P - p_c) * (460 + T_e))
            S_p = S_Th * (abs(base) ** 0.6)
            return [S_p]
        except:
            return []

    @staticmethod
    def eqn_10_19__T_e(
        P: float, S_Th: float, S_p: float, T_i: float, p_c: float, p_s: float
    ):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        try:
            power_term = (S_p / S_Th) ** (1/0.6)
            T_e = ((P - p_s) * (460 + T_i)) / ((P - p_c) * power_term) - 460
            return [T_e]
        except:
            return []

    @staticmethod
    def eqn_10_19__T_i(
        P: float, S_Th: float, S_p: float, T_e: float, p_c: float, p_s: float
    ):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        try:
            power_term = (S_p / S_Th) ** (1/0.6)
            T_i = (power_term * (P - p_c) * (460 + T_e)) / (P - p_s) - 460
            return [T_i]
        except:
            return []

    @staticmethod
    def eqn_10_19__p_c(
        P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_s: float
    ):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        try:
            power_term = (S_p / S_Th) ** (1/0.6)
            p_c = P - ((P - p_s) * (460 + T_i)) / ((460 + T_e) * power_term)
            return [p_c]
        except:
            return []

    @staticmethod
    def eqn_10_19__p_s(
        P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float
    ):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        try:
            power_term = (S_p / S_Th) ** (1/0.6)
            p_s = P - power_term * ((P - p_c) * (460 + T_e)) / (460 + T_i)
            return [p_s]
        except:
            return []


