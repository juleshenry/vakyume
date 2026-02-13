from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_6_11a(
        A_d: float = None,
        delta_T: float = None,
        delta_h_i: float = None,
        delta_m: float = None,
        h_d: float = None,
        m_b: float = None,
        t_R: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_6_11a__A_d(
        delta_T: float,
        delta_h_i: float,
        delta_m: float,
        h_d: float,
        m_b: float,
        t_R: float,
    ):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        A_d = delta_h_i * delta_m * m_b / (delta_T * h_d * t_R)
        result.append(A_d)
        return result

    @staticmethod
    def eqn_6_11a__delta_T(
        A_d: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float
    ):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_T = delta_h_i * delta_m * m_b / (A_d * h_d * t_R)
        result.append(delta_T)
        return result

    @staticmethod
    def eqn_6_11a__delta_h_i(
        A_d: float, delta_T: float, delta_m: float, h_d: float, m_b: float, t_R: float
    ):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_h_i = A_d * delta_T * h_d * t_R / (delta_m * m_b)
        result.append(delta_h_i)
        return result

    @staticmethod
    def eqn_6_11a__delta_m(
        A_d: float, delta_T: float, delta_h_i: float, h_d: float, m_b: float, t_R: float
    ):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_m = A_d * delta_T * h_d * t_R / (delta_h_i * m_b)
        result.append(delta_m)
        return result

    @staticmethod
    def eqn_6_11a__h_d(
        A_d: float,
        delta_T: float,
        delta_h_i: float,
        delta_m: float,
        m_b: float,
        t_R: float,
    ):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        h_d = delta_h_i * delta_m * m_b / (A_d * delta_T * t_R)
        result.append(h_d)
        return result

    @staticmethod
    def eqn_6_11a__m_b(
        A_d: float,
        delta_T: float,
        delta_h_i: float,
        delta_m: float,
        h_d: float,
        t_R: float,
    ):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        m_b = A_d * delta_T * h_d * t_R / (delta_h_i * delta_m)
        result.append(m_b)
        return result

    @staticmethod
    def eqn_6_11a__t_R(
        A_d: float,
        delta_T: float,
        delta_h_i: float,
        delta_m: float,
        h_d: float,
        m_b: float,
    ):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        t_R = delta_h_i * delta_m * m_b / (A_d * delta_T * h_d)
        result.append(t_R)
        return result




