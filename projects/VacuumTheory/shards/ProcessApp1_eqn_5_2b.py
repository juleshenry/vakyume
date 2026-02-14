from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_2b(K_1=None, K_2=None, x_1=None, x_2=None, y_1=None, y_2=None, **kwargs):
        return

    @staticmethod
    def eqn_5_2b__K_1(K_2: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_1 = K_2*x_2*y_1/(x_1*y_2)
        result.append(K_1)
        return result

    @staticmethod
    def eqn_5_2b__K_2(K_1: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_2 = K_1*x_1*y_2/(x_2*y_1)
        result.append(K_2)
        return result

    @staticmethod
    def eqn_5_2b__x_1(K_1: float, K_2: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_1 = K_2*x_2*y_1/(K_1*y_2)
        result.append(x_1)
        return result

    @staticmethod
    def eqn_5_2b__x_2(K_1: float, K_2: float, x_1: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_2 = K_1*x_1*y_2/(K_2*y_1)
        result.append(x_2)
        return result

    @staticmethod
    def eqn_5_2b__y_1(K_1: float, K_2: float, x_1: float, x_2: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_1 = K_1*x_1*y_2/(K_2*x_2)
        result.append(y_1)
        return result

    @staticmethod
    def eqn_5_2b__y_2(K_1: float, K_2: float, x_1: float, x_2: float, y_1: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_2 = K_2*x_2*y_1/(K_1*x_1)
        result.append(y_2)
        return result

