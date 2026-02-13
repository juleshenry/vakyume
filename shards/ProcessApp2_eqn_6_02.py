from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp2:
    @kwasak_static
    def eqn_6_02(
        Q_v: float = None,
        T_1: float = None,
        T_2: float = None,
        T_R: float = None,
        c_p: float = None,
        w_1: float = None,
        w_2: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_6_02__Q_v(
        T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        Q_v = c_p * (T_1 * w_1 + T_2 * w_2 - T_R * w_1 - T_R * w_2) / 12000
        result.append(Q_v)
        return result

    @staticmethod
    def eqn_6_02__T_1(
        Q_v: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_1 = (12000 * Q_v + T_R * c_p * w_1 + c_p * w_2 * (-T_2 + T_R)) / (c_p * w_1)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_02__T_2(
        Q_v: float, T_1: float, T_R: float, c_p: float, w_1: float, w_2: float
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_2 = (12000 * Q_v + T_R * c_p * w_2 + c_p * w_1 * (-T_1 + T_R)) / (c_p * w_2)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_02__T_R(
        Q_v: float, T_1: float, T_2: float, c_p: float, w_1: float, w_2: float
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_R = (-12000 * Q_v + T_1 * c_p * w_1 + T_2 * c_p * w_2) / (c_p * (w_1 + w_2))
        result.append(T_R)
        return result

    @staticmethod
    def eqn_6_02__c_p(
        Q_v: float, T_1: float, T_2: float, T_R: float, w_1: float, w_2: float
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        c_p = 12000 * Q_v / (T_1 * w_1 + T_2 * w_2 - T_R * w_1 - T_R * w_2)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_02__w_1(
        Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_2: float
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_1 = (12000 * Q_v - T_2 * c_p * w_2 + T_R * c_p * w_2) / (c_p * (T_1 - T_R))
        result.append(w_1)
        return result

    @staticmethod
    def eqn_6_02__w_2(
        Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_1: float
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_2 = (12000 * Q_v - T_1 * c_p * w_1 + T_R * c_p * w_1) / (c_p * (T_2 - T_R))
        result.append(w_2)
        return result


