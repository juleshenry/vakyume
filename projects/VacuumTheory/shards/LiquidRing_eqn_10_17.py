from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_17(P=None, S_0=None, S_Th=None, p_0=None, p_s=None, **kwargs):
        return

    @staticmethod
    def eqn_10_17__P(S_0: float, S_Th: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        P = (p_0*(S_Th/S_0)**(5/3) - p_s)/((S_Th/S_0)**1.66666666666667 - 1.0)
        result.append(P)
        P = (p_0*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - p_s)/((-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        P = (p_0*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - p_s)/((-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        return result

    @staticmethod
    def eqn_10_17__S_0(P: float, S_Th: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_0 = S_Th/((P - p_s)/(P - p_0))**(3/5)
        result.append(S_0)
        return result

    @staticmethod
    def eqn_10_17__S_Th(P: float, S_0: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_Th = S_0*((P - p_s)/(P - p_0))**(3/5)
        result.append(S_Th)
        return result

    @staticmethod
    def eqn_10_17__p_0(P: float, S_0: float, S_Th: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        p_0 = (P*(S_Th/S_0)**(5/3) - P + p_s)/(S_Th/S_0)**(5/3)
        result.append(p_0)
        p_0 = (P*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - P + p_s)/(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        p_0 = (P*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - P + p_s)/(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        return result

    @staticmethod
    def eqn_10_17__p_s(P: float, S_0: float, S_Th: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        p_s = -P*(S_Th/S_0)**(5/3) + P + p_0*(S_Th/S_0)**(5/3)
        result.append(p_s)
        p_s = -P*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 + P + p_0*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_s)
        p_s = -P*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 + P + p_0*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_s)
        return result

