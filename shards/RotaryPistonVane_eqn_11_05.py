from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class RotaryPistonVane:
    @kwasak_static
    def eqn_11_05(
        P_0_v: float = None,
        P_D: float = None,
        p_g: float = None,
        p_v_max: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_11_05__P_0_v(P_D: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_0_v = P_D * p_v_max / (p_g + p_v_max)
        result.append(P_0_v)
        return result

    @staticmethod
    def eqn_11_05__P_D(P_0_v: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_D = P_0_v * (p_g + p_v_max) / p_v_max
        result.append(P_D)
        return result

    @staticmethod
    def eqn_11_05__p_g(P_0_v: float, P_D: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_g = p_v_max * (-P_0_v + P_D) / P_0_v
        result.append(p_g)
        return result

    @staticmethod
    def eqn_11_05__p_v_max(P_0_v: float, P_D: float, p_g: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_v_max = -P_0_v * p_g / (P_0_v - P_D)
        result.append(p_v_max)
        return result


