from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_01(P: float = None, p_i: float = None, y_i: float = None, **kwargs):
        return

    @staticmethod
    def eqn_7_01__P(p_i: float, y_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        P = p_i / y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_7_01__p_i(P: float, y_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        p_i = P * y_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_01__y_i(P: float, p_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        y_i = p_i / P
        result.append(y_i)
        return result


