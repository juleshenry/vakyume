from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp2:
    @kwasak_static
    def eqn_6_01(
        T_1: float = None,
        T_2: float = None,
        T_R: float = None,
        c_p: float = None,
        del_h_v: float = None,
        w_1: float = None,
        w_2: float = None,
        w_v: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_6_01__T_1(
        T_2: float,
        T_R: float,
        c_p: float,
        del_h_v: float,
        w_1: float,
        w_2: float,
        w_v: float,
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_1 = (T_R * c_p * w_1 + c_p * w_2 * (-T_2 + T_R) + del_h_v * w_v) / (c_p * w_1)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_01__T_2(
        T_1: float,
        T_R: float,
        c_p: float,
        del_h_v: float,
        w_1: float,
        w_2: float,
        w_v: float,
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_2 = (T_R * c_p * w_2 + c_p * w_1 * (-T_1 + T_R) + del_h_v * w_v) / (c_p * w_2)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_01__T_R(
        T_1: float,
        T_2: float,
        c_p: float,
        del_h_v: float,
        w_1: float,
        w_2: float,
        w_v: float,
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_R = (T_1 * c_p * w_1 + T_2 * c_p * w_2 - del_h_v * w_v) / (c_p * (w_1 + w_2))
        result.append(T_R)
        return result

    @staticmethod
    def eqn_6_01__c_p(
        T_1: float,
        T_2: float,
        T_R: float,
        del_h_v: float,
        w_1: float,
        w_2: float,
        w_v: float,
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        c_p = del_h_v * w_v / (T_1 * w_1 + T_2 * w_2 - T_R * w_1 - T_R * w_2)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_01__del_h_v(
        T_1: float,
        T_2: float,
        T_R: float,
        c_p: float,
        w_1: float,
        w_2: float,
        w_v: float,
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        del_h_v = c_p * (T_1 * w_1 + T_2 * w_2 - T_R * w_1 - T_R * w_2) / w_v
        result.append(del_h_v)
        return result

    @staticmethod
    def eqn_6_01__w_1(
        T_1: float,
        T_2: float,
        T_R: float,
        c_p: float,
        del_h_v: float,
        w_2: float,
        w_v: float,
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_1 = (-T_2 * c_p * w_2 + T_R * c_p * w_2 + del_h_v * w_v) / (c_p * (T_1 - T_R))
        result.append(w_1)
        return result

    @staticmethod
    def eqn_6_01__w_2(
        T_1: float,
        T_2: float,
        T_R: float,
        c_p: float,
        del_h_v: float,
        w_1: float,
        w_v: float,
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_2 = (-T_1 * c_p * w_1 + T_R * c_p * w_1 + del_h_v * w_v) / (c_p * (T_2 - T_R))
        result.append(w_2)
        return result

    @staticmethod
    def eqn_6_01__w_v(
        T_1: float,
        T_2: float,
        T_R: float,
        c_p: float,
        del_h_v: float,
        w_1: float,
        w_2: float,
    ):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_v = c_p * (T_1 * w_1 + T_2 * w_2 - T_R * w_1 - T_R * w_2) / del_h_v
        result.append(w_v)
        return result


