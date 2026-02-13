from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_4_10(
        T: float = None,
        V: float = None,
        del_P: float = None,
        leakage: float = None,
        t: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_4_10__T(V: float, del_P: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        T = 3.127 * V * del_P / (leakage * t)
        result.append(T)
        return result

    @staticmethod
    def eqn_4_10__V(T: float, del_P: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        V = 0.319795330988168 * T * leakage * t / del_P
        result.append(V)
        return result

    @staticmethod
    def eqn_4_10__del_P(T: float, V: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        del_P = 0.319795330988168 * T * leakage * t / V
        result.append(del_P)
        return result

    @staticmethod
    def eqn_4_10__leakage(T: float, V: float, del_P: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        leakage = 3.127 * V * del_P / (T * t)
        result.append(leakage)
        return result

    @staticmethod
    def eqn_4_10__t(T: float, V: float, del_P: float, leakage: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        t = 3.127 * V * del_P / (T * leakage)
        result.append(t)
        return result




