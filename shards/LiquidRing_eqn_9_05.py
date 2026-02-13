from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_9_05(
        V: float = None,
        r_h: float = None,
        t_h: float = None,
        w_h: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_9_05__V(r_h: float, t_h: float, w_h: float):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        V = t_h * w_h / r_h
        result.append(V)
        return result

    @staticmethod
    def eqn_9_05__r_h(V: float, t_h: float, w_h: float):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        r_h = t_h * w_h / V
        result.append(r_h)
        return result

    @staticmethod
    def eqn_9_05__t_h(V: float, r_h: float, w_h: float):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        t_h = V * r_h / w_h
        result.append(t_h)
        return result

    @staticmethod
    def eqn_9_05__w_h(V: float, r_h: float, t_h: float):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        w_h = V * r_h / t_h
        result.append(w_h)
        return result




