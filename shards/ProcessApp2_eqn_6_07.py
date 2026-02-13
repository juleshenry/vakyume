from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp2:
    @kwasak_static
    def eqn_6_07(
        C_1: float = None,
        C_2: float = None,
        T_1: float = None,
        T_2: float = None,
        c_p: float = None,
        delta_h_c: float = None,
        delta_h_v: float = None,
        m_b: float = None,
        m_v: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_6_07__C_1(
        C_2: float,
        T_1: float,
        T_2: float,
        c_p: float,
        delta_h_c: float,
        delta_h_v: float,
        m_b: float,
        m_v: float,
    ):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_1 = (C_2 * delta_h_c * m_b + c_p * m_b * (-T_1 + T_2) + delta_h_v * m_v) / (
            delta_h_c * m_b
        )
        result.append(C_1)
        return result

    @staticmethod
    def eqn_6_07__C_2(
        C_1: float,
        T_1: float,
        T_2: float,
        c_p: float,
        delta_h_c: float,
        delta_h_v: float,
        m_b: float,
        m_v: float,
    ):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_2 = (C_1 * delta_h_c * m_b + c_p * m_b * (T_1 - T_2) - delta_h_v * m_v) / (
            delta_h_c * m_b
        )
        result.append(C_2)
        return result

    @staticmethod
    def eqn_6_07__T_1(
        C_1: float,
        C_2: float,
        T_2: float,
        c_p: float,
        delta_h_c: float,
        delta_h_v: float,
        m_b: float,
        m_v: float,
    ):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_1 = (T_2 * c_p * m_b + delta_h_c * m_b * (-C_1 + C_2) + delta_h_v * m_v) / (
            c_p * m_b
        )
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_07__T_2(
        C_1: float,
        C_2: float,
        T_1: float,
        c_p: float,
        delta_h_c: float,
        delta_h_v: float,
        m_b: float,
        m_v: float,
    ):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_2 = (T_1 * c_p * m_b + delta_h_c * m_b * (C_1 - C_2) - delta_h_v * m_v) / (
            c_p * m_b
        )
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_07__c_p(
        C_1: float,
        C_2: float,
        T_1: float,
        T_2: float,
        delta_h_c: float,
        delta_h_v: float,
        m_b: float,
        m_v: float,
    ):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        c_p = (-C_1 * delta_h_c * m_b + C_2 * delta_h_c * m_b + delta_h_v * m_v) / (
            m_b * (T_1 - T_2)
        )
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_07__delta_h_c(
        C_1: float,
        C_2: float,
        T_1: float,
        T_2: float,
        c_p: float,
        delta_h_v: float,
        m_b: float,
        m_v: float,
    ):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_c = (-T_1 * c_p * m_b + T_2 * c_p * m_b + delta_h_v * m_v) / (
            m_b * (C_1 - C_2)
        )
        result.append(delta_h_c)
        return result

    @staticmethod
    def eqn_6_07__delta_h_v(
        C_1: float,
        C_2: float,
        T_1: float,
        T_2: float,
        c_p: float,
        delta_h_c: float,
        m_b: float,
        m_v: float,
    ):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_v = (
            m_b * (C_1 * delta_h_c - C_2 * delta_h_c + T_1 * c_p - T_2 * c_p) / m_v
        )
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_07__m_b(
        C_1: float,
        C_2: float,
        T_1: float,
        T_2: float,
        c_p: float,
        delta_h_c: float,
        delta_h_v: float,
        m_v: float,
    ):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_b = (
            delta_h_v
            * m_v
            / (C_1 * delta_h_c - C_2 * delta_h_c + T_1 * c_p - T_2 * c_p)
        )
        result.append(m_b)
        return result

    @staticmethod
    def eqn_6_07__m_v(
        C_1: float,
        C_2: float,
        T_1: float,
        T_2: float,
        c_p: float,
        delta_h_c: float,
        delta_h_v: float,
        m_b: float,
    ):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_v = (
            m_b
            * (C_1 * delta_h_c - C_2 * delta_h_c + T_1 * c_p - T_2 * c_p)
            / delta_h_v
        )
        result.append(m_v)
        return result


