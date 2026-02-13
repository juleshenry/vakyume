from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_01(K_i: float = None, x_i: float = None, y_i: float = None, **kwargs):
        return

    @staticmethod
    def eqn_5_01__K_i(x_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        K_i = y_i / x_i
        result.append(K_i)
        return result

    @staticmethod
    def eqn_5_01__x_i(K_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        x_i = y_i / K_i
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_01__y_i(K_i: float, x_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        y_i = K_i * x_i
        result.append(y_i)
        return result


