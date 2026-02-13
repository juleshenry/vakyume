from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_03(
        P_i_0: float = None,
        epsilon_i: float = None,
        p_i: float = None,
        x_i: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_7_03__P_i_0(epsilon_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        P_i_0 = p_i / (epsilon_i * x_i)
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_03__epsilon_i(P_i_0: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        epsilon_i = p_i / (P_i_0 * x_i)
        result.append(epsilon_i)
        return result

    @staticmethod
    def eqn_7_03__p_i(P_i_0: float, epsilon_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        p_i = P_i_0 * epsilon_i * x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_03__x_i(P_i_0: float, epsilon_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        x_i = p_i / (P_i_0 * epsilon_i)
        result.append(x_i)
        return result


