from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_15(P=None, S_Th=None, S_p=None, p_s=None, **kwargs):
        return

    @staticmethod
    def eqn_10_15__P(S_Th: float, S_p: float, p_s: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        P = S_Th*p_s/(S_Th - S_p)
        result.append(P)
        return result

    @staticmethod
    def eqn_10_15__S_Th(P: float, S_p: float, p_s: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        S_Th = P*S_p/(P - p_s)
        result.append(S_Th)
        return result

    @staticmethod
    def eqn_10_15__S_p(P: float, S_Th: float, p_s: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        S_p = S_Th*(P - p_s)/P
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_15__p_s(P: float, S_Th: float, S_p: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        p_s = P*(S_Th - S_p)/S_Th
        result.append(p_s)
        return result

